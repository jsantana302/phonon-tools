#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# HIPHIVE wrapper for CP2K
################################################################################
# Copyright:
#   Julia Santana (2025)
#
# License:
#   GNU General Public License v3 or later
#   <http://www.gnu.org/licenses/>
################################################################################

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from ase.io import read, write
from ase import Atoms
from hiphive.cutoffs import BaseClusterFilter
from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
from hiphive import ForceConstants
from trainstation import CrossValidationEstimator

# -----------------------------------------------------------------------------#
# Globals
# -----------------------------------------------------------------------------#

BINS = {
    "displacement": np.linspace(0.0, 1.5, 300),
    "distance": np.linspace(1.0, 5.0, 300),
}

HA_PER_BOHR_TO_EV_PER_ANG = 27.211386245988 / 0.529177210903  # ≈ 51.422067


# -----------------------------------------------------------------------------#
# Basic helpers
# -----------------------------------------------------------------------------#

def fix_format(seq: Sequence) -> str:
    """Render a Python sequence as space-separated list (CONTROL cards)."""
    return " ".join(str(x) for x in seq)


def get_histogram_data(
    data: Iterable[float],
    bins: int | Sequence[float] = 100,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return normalized histogram centers and densities."""
    arr = list(data)
    counts, edges = np.histogram(arr, bins=bins, density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    return centers, counts


def get_distributions(structure_list: List[Atoms], ref_pos: np.ndarray
                      ) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Compute distributions of interatomic distances and displacements."""
    distances, displacements = [], []
    for atoms in structure_list:
        distances.extend(atoms.get_all_distances(mic=True).ravel())
        displacements.extend(np.linalg.norm(atoms.positions - ref_pos, axis=1))
    print("Max displacement in distorted structures:",
          f"{np.max(displacements):.4f} Å")
    return {
        "distance": get_histogram_data(distances, BINS["distance"]),
        "displacement": get_histogram_data(displacements, BINS["displacement"]),
    }


def save_json(obj, path: str) -> None:
    """Write a small JSON file with UTF-8 encoding."""
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


# -----------------------------------------------------------------------------#
# CP2K I/O helpers
# -----------------------------------------------------------------------------#

def cp2k_dir_finished(directory: Path) -> bool:
    """Return True if a CP2K .out indicates a normal end.

    Heuristics: look for lines like "PROGRAM ENDED AT" or "CP2K ended
    normally" in any *.out in the directory.
    """
    for out in sorted(directory.glob("*.out")):
        try:
            t = out.read_text(errors="ignore")
        except Exception:
            continue
        if ("PROGRAM ENDED AT" in t) or ("CP2K ended normally" in t):
            return True
    return False





def write_cp2k_dirs(structures: List[Atoms],
                    template_dir: str,
                    coord_filename: str = "coord.xyz",
                    inp_glob: str = "*.inp") -> None:
    """Create job folders 00, 01, 02, ... in CWD, copy template, write coord.xyz, patch inputs.

    New layout (preferred): ./00, ./01, ...
    Legacy layout still read by read_cp2k_outputs: CALC/job-XX
    """
    template_dir = Path(template_dir)
    if not template_dir.is_dir():
        raise FileNotFoundError(f"Template dir not found: {template_dir}")

    # Backup any existing two-digit numeric job folders
    existing = [p for p in Path(".").iterdir()
                if p.is_dir() and re.fullmatch(r"\d{2}", p.name)]
    if existing:
        ts = datetime.now().strftime("%d-%b-%Y(%H:%M.%S)")
        backup_root = Path(f"JOBS_OLD_{ts}")
        print(f"Moving existing two-digit job folders -> {backup_root}/")
        backup_root.mkdir()
        for p in sorted(existing, key=lambda x: int(x.name)):
            shutil.move(str(p), str(backup_root / p.name))

    for i, structure in enumerate(structures, start=0):
        job_dir = Path(f"{i:02d}")
        shutil.copytree(template_dir, job_dir)
        # coord.xyz
        write(job_dir / coord_filename, structure, format="xyz")
        print(f"Created {job_dir}")


def _read_forces_from_forces_xyz(path: Path) -> np.ndarray | None:
    try:
        ats = read(path, index=-1)  # last frame
    except Exception:
        return None
    if "forces" in ats.arrays:
        return np.asarray(ats.arrays["forces"], float)
    fx, fy, fz = ats.arrays.get("fx"), ats.arrays.get("fy"), ats.arrays.get("fz")
    if fx is not None and fy is not None and fz is not None:
        return np.vstack([fx, fy, fz]).T.astype(float)
    return None

def _parse_forces_from_cp2k_out(out_path: Path, natoms: int) -> np.ndarray | None:
    """Parse 'ATOMIC FORCES' table from CP2K .out and return eV/Å."""
    try:
        txt = out_path.read_text(errors="ignore")
    except Exception:
        return None

    # Locate "ATOMIC FORCES in [units]"
    m = re.search(r"^\s*ATOMIC FORCES.*\[(?P<u>.+?)\].*$", txt, re.M)
    if not m:
        return None
    units = m.group("u").strip().lower()  # e.g., 'a.u.' or 'ev/angstrom'

    # Find block start line index
    lines = txt.splitlines()
    start_idx = None
    for idx, L in enumerate(lines):
        if re.match(r"^\s*ATOMIC FORCES", L):
            start_idx = idx
            break
    if start_idx is None:
        return None

    # Usually there are one or two header lines; parse natoms lines after
    # from the first line that looks like: index  element  fx   fy   fz
    data: List[List[float]] = []
    for L in lines[start_idx + 1: start_idx + 1 + natoms + 3]:
        if not L.strip():
            continue
        # extract last three floats on the line
        floats = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", L)
        if len(floats) >= 3:
            try:
                fx, fy, fz = map(float, floats[-3:])
            except Exception:
                continue
            data.append([fx, fy, fz])
            if len(data) == natoms:
                break

    if len(data) != natoms:
        return None

    F = np.array(data, float)
    if "a.u" in units or "hartree/bohr" in units:
        F = F * HA_PER_BOHR_TO_EV_PER_ANG
    # else assume already eV/Å
    return F


def read_cp2k_outputs(sc: Atoms,
                      rat_n: int,
                      coord_filename: str = "coord.xyz") -> List[Atoms]:
    """Collect forces/displacements from CP2K jobs into a structure list.

    Strategy
    --------
    Preferred new layout:
      - Iterate ./00, ./01, ... (up to `rat_n`).
    Legacy fallback:
      - Iterate CALC/job-XX (up to `rat_n`).

    Positions: read `coord.xyz` written by the wrapper.
    Forces: prefer `forces.xyz`; fallback to parsing any *.out.
    Build a copy of `sc` with arrays: 'displacements', 'forces'.

    Returns
    -------
    structures : list of ase.Atoms
    """
    from hiphive.utilities import get_displacements

    # Discover job directories (new style preferred)
    two_digit_dirs = sorted(
        [p for p in Path(".").iterdir() if p.is_dir() and re.fullmatch(r"\d{2}", p.name)],
        key=lambda p: int(p.name)
    )[:rat_n]

    job_dirs: List[Path]
    if two_digit_dirs:
        job_dirs = two_digit_dirs
    else:
        base = Path("CALC")
        if not base.is_dir():
            raise FileNotFoundError(
                "No two-digit job folders found and CALC/ does not exist. "
                "Run -m pre first?"
            )
        # legacy: CALC/job-XX (sort by trailing number if present)
        legacy = [p for p in base.glob("job-*") if p.is_dir()]
        def _legacy_key(p: Path):
            m = re.search(r"(\d+)$", p.name)
            return int(m.group(1)) if m else 10**9
        job_dirs = sorted(legacy, key=_legacy_key)[:rat_n]
        if not job_dirs:
            raise RuntimeError("No job-* directories found in CALC/.")

    structures: List[Atoms] = []
    all_force_components: List[float] = []

    for i, job in enumerate(job_dirs, start=1):
        # positions (the distorted structure used)
        cpath = job / coord_filename
        if not cpath.is_file():
            raise FileNotFoundError(f"Missing {cpath} (template must use it).")
        displaced = read(cpath)
        # ensure same cell and PBC as sc
        at = sc.copy()
        at.set_positions(displaced.get_positions())

        # forces: try forces.xyz, then parse .out
        F = None
        cand_paths = []
        for pat in ["forces.xyz", "force.xyz", "cp2k_forces.xyz",
                    "*forces*.xyz", "*FORCES*.xyz", "*Force*.xyz"]:
            cand_paths.extend(sorted(job.glob(pat)))

        for fpath in cand_paths:
            # try Extended XYZ
            F = _read_forces_from_forces_xyz(fpath)
            if F is None:
                # NEW: fall back to CP2K table parser on this same file
                F = _parse_forces_from_cp2k_out(fpath, natoms=len(sc))
            if F is not None:
                break
            
        if F is None:
            # fallback: parse any *.out
            for op in sorted(job.glob("*.out")):
                F = _parse_forces_from_cp2k_out(op, natoms=len(sc))
                if F is not None:
                    break
        if F is None:
            raise RuntimeError(f"Could not find forces in {job}")

                # displacements via tool (handles MIC)
        D = get_displacements(at, sc)

        # sanity
        if np.linalg.norm(D, axis=1).max() >= 1.0:
            raise AssertionError("Displacements >= 1.0 Å detected.")

        at_tmp = sc.copy()
        at_tmp.new_array("displacements", D)
        at_tmp.new_array("forces", F)
        structures.append(at_tmp)
        all_force_components.extend(np.abs(F).ravel())

        if not cp2k_dir_finished(job):
            print(f"[warn] CP2K did not confirm normal end in {job}")

        print(f"[{i:02d}] max|F|={np.max(np.abs(F)):.6f} eV/Å")

    if all_force_components:
        print(f"Global max |F|: {np.max(all_force_components):.6f} eV/Å")
        print(f"Global mean |F|: {np.mean(all_force_components):.6f} eV/Å")

    # export for provenance
    out = Path("rattled_structures.extxyz")
    if out.exists():
        out.unlink()
    write(out, structures)
    return structures


# -----------------------------------------------------------------------------#
# Hiphive-specific helpers
# -----------------------------------------------------------------------------#

def get_index_offset(atoms, atoms_ref) -> Tuple[List[int], List[np.ndarray]]:
    """Map indices from a repeated structure to the primitive reference."""
    from hiphive.core.atoms import spos_to_atom
    basis = atoms_ref.get_scaled_positions()
    spos = atoms_ref.cell.scaled_positions(atoms.positions)
    indices, offsets = [], []
    for s in spos:
        a = spos_to_atom(s, basis)
        indices.append(a.site)
        offsets.append(a.offset)
    return indices, offsets


class MaxCut2nd(BaseClusterFilter):
    """Cluster filter limiting 2nd-order pair vectors component-wise."""

    def __init__(self, lim: Sequence[float]):
        self.limx, self.limy, self.limz = lim

    def setup(self, atoms):
        self.atoms = atoms

    def __call__(self, cluster) -> bool:
        if len(cluster) >= 3:
            return True
        i, j = cluster[0], cluster[1]
        vec = self.atoms.get_distance(i, j, mic=True, vector=True)
        return (abs(vec[0]) < self.limx and
                abs(vec[1]) < self.limy and
                abs(vec[2]) < self.limz)


def collect_force_components(structures) -> np.ndarray:
    """Flatten all |F| components from a list of structures."""
    arrs = []
    for s in structures:
        if "forces" not in s.arrays:
            raise ValueError("Structure missing 'forces' array.")
        arrs.append(np.abs(np.asanyarray(s.arrays["forces"])).ravel())
    return np.concatenate(arrs) if arrs else np.array([])


def trim_outlier_structures(structures, k: float) -> list:
    """Drop structures with max |F| component > (Q3 + k·IQR)."""
    if k <= 0.0:
        return structures
    comps = collect_force_components(structures)
    if comps.size == 0:
        return structures
    q1, q3 = np.quantile(comps, [0.25, 0.75])
    iqr = q3 - q1
    thresh = q3 + k * iqr
    kept, dropped = [], 0
    for s in structures:
        max_comp = float(np.max(np.abs(s.arrays["forces"])))
        if max_comp <= thresh:
            kept.append(s)
        else:
            dropped += 1
    print(f"Outlier trim (k={k:.3f}): threshold={thresh:.6f} eV/Å, "
          f"dropped={dropped}, kept={len(kept)}")
    return kept


def build_cutoffs(sc: Atoms, cutoffs_raw: Sequence[str | float]
                  ) -> Tuple[List[float], bool, List[float]]:
    """Build numeric cutoffs from raw user entries."""
    folding = False
    lim: List[float] = []

    tok = str(cutoffs_raw[0]).upper() if cutoffs_raw else ""
    if tok == "MAX":
        border = [np.linalg.norm(sc.cell[i] / 2.0) - 1e-2 for i in range(3)]
        cutoffs_raw = [min(border)] + list(cutoffs_raw[1:])
    elif tok == "MAX2":
        lim = [np.linalg.norm(sc.cell[i] / 2.0) - 1e-2 for i in range(3)]
        c0 = math.sqrt(sum(x * x for x in lim))
        cutoffs_raw = [c0] + list(cutoffs_raw[1:])
    elif tok == "MAX3":
        lim = [np.linalg.norm(sc.cell[i] / 2.0) + 1e-2 for i in range(3)]
        c0 = math.sqrt(sum(x * x for x in lim))
        cutoffs_raw = [c0] + list(cutoffs_raw[1:])
        folding = True

    cutoffs: List[float] = []
    for x in cutoffs_raw:
        val = float(x)
        if val > 0.0:
            cutoffs.append(val)
        else:
            from hiphive.utilities import get_neighbor_shells
            shells = get_neighbor_shells(sc, 15)
            index = int(abs(val)) + 1
            value = shells[index].distance - 0.05
            cutoffs.append(value)
    return cutoffs, folding, lim


# -----------------------------------------------------------------------------#
# Mode runners
# -----------------------------------------------------------------------------#

def run_control(ngrid: Sequence[str], scell: Sequence[str], temp: float) -> None:
    """Write ShengBTE CONTROL from 'POSCAR'/'pc.xyz' in cwd (generic)."""
    atoms = read("POSCAR") if Path("POSCAR").exists() else read("pc.xyz")
    all_symbols = atoms.get_chemical_symbols()
    elements = list(dict.fromkeys(all_symbols))

    with open("CONTROL", "w", encoding="utf-8") as f:
        f.write("&allocations\n")
        f.write(f"\t nelements={len(elements)},\n")
        f.write(f"\t natoms={len(all_symbols)},\n")
        f.write(f"\t ngrid(:)={' '.join(ngrid)}\n")
        f.write("&end\n")

        f.write("&crystal\n\t lfactor=0.1,\n")
        for i, cell in enumerate(atoms.cell, start=1):
            f.write(
                f"\t lattvec(:,{i})="
                f"{cell[0]:.12f}\t{cell[1]:.12f}\t{cell[2]:.12f},\n"
            )
        f.write(f"\t elements={fix_format(elements)},\n")

        idx = {el: i + 1 for i, el in enumerate(elements)}
        types = [idx[s] for s in all_symbols]
        f.write(f"\t types={fix_format(types)},\n")

        for i, s in enumerate(atoms.get_scaled_positions(), start=1):
            f.write(
                f"\t positions(:,{i})="
                f"{s[0]:.12f}\t{s[1]:.12f}\t{s[2]:.12f},\n"
            )
        f.write(f"\t scell(:)={' '.join(scell)}\n&end\n")

        f.write("&parameters\n")
        f.write(f"\tT={temp}\n\tscalebroad=0.1\n\tmaxiter=50\n&end\n")
        f.write("&flags\n\tisotopes=.TRUE.\n\tnonanalytic=.FALSE.\n"
                "\tnanowires=.FALSE.\n&end\n")


def run_pre(sc: Atoms, rat_mode: str, rat_n: int, rat_std: float,
            template_dir: str, coord_filename: str, inp_glob: str) -> None:
    """Generate distorted structures and CP2K folders."""
    if rat_mode == "rattle":
        from hiphive.structure_generation import generate_rattled_structures
        from numpy.random import randint
        seed = randint(1)
        structures = generate_rattled_structures(sc, rat_n, rat_std, seed=seed)

    elif rat_mode == "mc":
        from hiphive.structure_generation import generate_mc_rattled_structures
        from numpy.random import randint
        from hiphive.utilities import get_neighbor_shells
        seed = randint(1)
        shells = get_neighbor_shells(sc, 6)
        dmin = shells[0].distance - 0.1
        structures = generate_mc_rattled_structures(
            sc, rat_n, rat_std * 0.15, dmin, seed=seed, n_iter=20
        )

    elif rat_mode == "phon-rat":
        from hiphive.structure_generation import generate_phonon_rattled_structures
        from hiphive import ForceConstantPotential
        temp = 300.0
        fcp = ForceConstantPotential.read("Potential.fcp")
        fc2 = fcp.get_force_constants(sc).get_fc_array(order=2, format="ase")
        structures = generate_phonon_rattled_structures(
            sc, fc2, rat_n, temp, imag_freq_factor=2
        )
    else:
        raise ValueError(f"Unknown rat_mode: {rat_mode!r}")

    write_cp2k_dirs(structures, template_dir, coord_filename, inp_glob)

    ref_pos = sc.get_positions()
    distributions = get_distributions(structures, ref_pos)
    disp_x, disp_y = distributions["displacement"]
    dist_x, dist_y = distributions["distance"]

    with open("histogram.dat", "w", encoding="utf-8") as f:
        f.write(f"{'Ampl':>10s}  {'Histo':>10s}  {'r':>10s}  {'g(r)':>10s}\n")
        for a, b, c, d in zip(disp_x, disp_y, dist_x, dist_y):
            f.write(f"{a:10.4f}  {b:10.4f}  {c:10.4f}  {d:10.4f}\n")

    if any(a > 0.8 and b > 0.0 for a, b in zip(disp_x, disp_y)):
        best = float(disp_x[int(np.argmax(disp_y))])
        print("WARNING: Displacements > 0.8 Å detected; consider MC rattle "
              f"with amplitude ≈ {best:.3f} Å if appropriate.")


def run_shell(sc: Atoms) -> None:
    """Print neighbor shells (debug helper)."""
    from hiphive.utilities import get_neighbor_shells
    for i, sh in enumerate(get_neighbor_shells(sc, 15)):
        print(i, sh)


def run_post(pc: Atoms, sc: Atoms, cutoffs: List[float], folding: bool,
             lim: List[float], rat_n: int, fit_methods: Sequence[str],
             train_fraction: float, benchmark: bool, rot_sumrule: bool,
             trim_iqr_k: float = 0.0,
             cutoff_sweep: Sequence[float] | None = None,
             bootstrap: int = 0,
             report_json: str = "fit_report.json",
             coord_filename: str = "coord.xyz") -> None:
    """Fit FCs and export files, with trimming/sweep/bootstrap/report."""


    print(f"Reading up to {rat_n} CP2K jobs …")
    structures_all = read_cp2k_outputs(sc, rat_n, coord_filename=coord_filename)
    structures_used = trim_outlier_structures(structures_all, trim_iqr_k)

    def _train_one(current_cutoffs, tag_prefix: str, frames):
        results = []
        print(f"Cutoffs = {current_cutoffs}")

        if lim:
            cs = ClusterSpace(sc, current_cutoffs, cluster_filter=MaxCut2nd(lim))
        else:
            cs = ClusterSpace(sc, current_cutoffs)

        stc = StructureContainer(cs)
        if folding:
            sc2 = sc.repeat(2)
            for s in frames:
                stc.add_structure(s.repeat(2))
        else:
            for s in frames:
                stc.add_structure(s)

        sizes = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8] if benchmark \
            else [train_fraction]

        fit_kwargs = defaultdict(dict)
        fit_kwargs["lasso"] = dict(max_iter=5000)

        for size in sizes:
            for method in fit_methods:
                print(f"[fit] method={method} train_size={size}")
                cve = CrossValidationEstimator(
                    stc.get_fit_data(),
                    fit_method=method,
                    validation_method="shuffle-split",
                    train_size=size,
                    test_size=1 - size,
                    n_splits=5,
                    **fit_kwargs.get(method, {}),
                )
                cve.validate()
                base_splits = list(cve.rmse_validation_splits)
                cve.train()

                # extra bootstrap validations
                boot_splits = []
                if bootstrap > 0:
                    cve_boot = CrossValidationEstimator(
                        stc.get_fit_data(),
                        fit_method=method,
                        validation_method="shuffle-split",
                        train_size=size,
                        test_size=1 - size,
                        n_splits=bootstrap,
                        **fit_kwargs.get(method, {}),
                    )
                    cve_boot.validate()
                    boot_splits = list(cve_boot.rmse_validation_splits)

                # build FCP
                if rot_sumrule:
                    from hiphive import enforce_rotational_sum_rules
                    pars = enforce_rotational_sum_rules(
                        cs, cve.parameters, ["Huang", "Born-Huang"]
                    )
                else:
                    pars = cve.parameters
                fcp = ForceConstantPotential(cs, pars)

                label = "_" + "_".join(f"{x:.4f}" for x in current_cutoffs) + "_"
                dirname = f"{tag_prefix}cutoff{label}n_{rat_n}_{method}_tf_{size}"
                Path(dirname).mkdir(exist_ok=False)
                fcp.write("Potential.fcp")
                os.replace("Potential.fcp", Path(dirname) / "Potential.fcp")

                # FCs
                if folding:
                    fcs2 = fcp.get_force_constants(sc2)
                    fc2 = fcs2.get_fc_array(order=2)
                    fcs2_2 = ForceConstants.from_arrays(sc2, fc2_array=fc2,
                                                        fc3_array=None)
                    n_atoms = len(sc)
                    indices, _ = get_index_offset(sc2, sc)
                    fc2_fold = np.zeros((n_atoms, n_atoms, 3, 3))
                    for i in range(n_atoms):
                        for i2 in range(len(sc2)):
                            j = indices[i2]
                            fc2_fold[i, j] += fcs2_2[i, i2]
                    fcs = ForceConstants.from_arrays(sc, fc2_array=fc2_fold,
                                                     fc3_array=None)
                else:
                    fcs = fcp.get_force_constants(sc)

                # export
                fcs.write_to_phonopy("FORCE_CONSTANTS_2ND", format="text")
                fcs.write_to_shengBTE("FORCE_CONSTANTS_3RD", pc)
                os.replace("FORCE_CONSTANTS_2ND",
                           Path(dirname) / "FORCE_CONSTANTS_2ND")
                os.replace("FORCE_CONSTANTS_3RD",
                           Path(dirname) / "FORCE_CONSTANTS_3RD")

                entry = {
                    "dir": dirname,
                    "method": method,
                    "train_size": float(size),
                    "cutoffs": list(map(float, current_cutoffs)),
                    "cv_rmse_splits": base_splits,
                    "cv_rmse_mean": float(np.mean(base_splits)),
                    "cv_rmse_std": float(np.std(base_splits)),
                    "bootstrap_splits": boot_splits,
                    "bootstrap_mean": (
                        float(np.mean(boot_splits)) if boot_splits else None
                    ),
                    "bootstrap_std": (
                        float(np.std(boot_splits)) if boot_splits else None
                    ),
                    "rot_sumrule": bool(rot_sumrule),
                    "folding": bool(folding),
                    "n_structures": int(len(frames)),
                }
                save_json(entry, str(Path(dirname) / "fit_report.json"))
                results.append(entry)
        return results

    all_results = {
        "sweep": [],
        "base_cutoffs": list(map(float, cutoffs)),
        "trim_iqr_k": float(trim_iqr_k),
        "bootstrap": int(bootstrap),
        "benchmark": bool(benchmark),
    }

    if cutoff_sweep:
        print("Cutoff sweep over first cutoff:", cutoff_sweep)
        for c0 in cutoff_sweep:
            c0 = float(c0)
            current = [c0] + list(cutoffs[1:])
            tag = f"DIR_sweep_c0_{c0:.3f}_"
            res = _train_one(current, tag, structures_used)
            all_results["sweep"].append({"c0": c0, "results": res})
    else:
        res = _train_one(cutoffs, "DIR_", structures_used)
        all_results["sweep"].append({"c0": float(cutoffs[0]), "results": res})

    save_json(all_results, report_json)
    print(f"Wrote aggregate report: {report_json}")


def run_scph(sc: Atoms, cutoffs: List[float], temps: Sequence[float],
             alpha: float, n_iterations: int, n_structures: int,
             rot_sumrule: bool) -> None:
    """Self-consistent harmonic loop and export FC2(T)."""
    from hiphive.self_consistent_phonons import \
        self_consistent_harmonic_model
    from hiphive.calculators import ForceConstantCalculator
    from hiphive import ClusterSpace, ForceConstantPotential

    print("Reading Potential.fcp …")
    fcp = ForceConstantPotential.read("Potential.fcp")
    cs2 = ClusterSpace(sc, [cutoffs[0]])
    fcs0 = fcp.get_force_constants(sc)
    calc = ForceConstantCalculator(fcs0)

    Path("scph_trajs").mkdir(exist_ok=True)
    parameters_start = None
    for T in temps:
        print(f"Temperature = {T} K")
        parameters_traj = self_consistent_harmonic_model(
            sc, calc, cs2, T, alpha, n_iterations, n_structures,
            imag_freq_factor=2, parameters_start=parameters_start
        )
        parameters_start = parameters_traj[-1]

        if rot_sumrule:
            from hiphive import enforce_rotational_sum_rules
            pars_rot = enforce_rotational_sum_rules(
                cs2, parameters_traj[-1], ["Huang", "Born-Huang"]
            )
            fcp_scph = ForceConstantPotential(cs2, pars_rot)
        else:
            fcp_scph = ForceConstantPotential(cs2, parameters_traj[-1])

        fcp_scph.write(f"Potential_T{T}.fcp")
        np.savetxt(f"scph_trajs/scph_parameters_T{T}",
                   np.array(parameters_traj))

        fcs_scph = fcp_scph.get_force_constants(sc)
        fcs_scph.write_to_phonopy(f"FORCE_CONSTANTS_2ND_scph_{T}",
                                  format="text")


def run_ehm(sc: Atoms, temps: Sequence[float], size: int,
            neq: int, nprod: int, dt_fs: float, dump: int,
            cutoffs: List[float]) -> None:
    """Effective harmonic model via MD with FC calculator."""
    from ase import units
    from ase.io.trajectory import Trajectory
    from ase.md import MDLogger
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from hiphive import ClusterSpace, ForceConstantPotential
    from hiphive.calculators import ForceConstantCalculator
    from hiphive import StructureContainer
    from trainstation import Optimizer

    Path("md_runs").mkdir(exist_ok=True)

    print("Reading Potential.fcp …")
    fcp = ForceConstantPotential.read("Potential.fcp")
    prim = fcp.primitive_structure

    atoms_ideal = prim.repeat(size)
    fcs = fcp.get_force_constants(atoms_ideal)
    calc = ForceConstantCalculator(fcs)

    for T in temps:
        print(f"Temperature: {T} K")

        atoms = atoms_ideal.copy()
        atoms.set_calculator(calc)
        dyn = Langevin(atoms, dt_fs * units.fs, T * units.kB, 0.02)

        rs = np.random.RandomState(2020)
        MaxwellBoltzmannDistribution(atoms, T * units.kB, rng=rs)
        dyn.run(neq)

        log_file = f"md_runs/T{T}.log"
        traj_file = f"md_runs/T{T}.traj"
        logger = MDLogger(dyn, atoms, log_file, header=True,
                          stress=False, peratom=True, mode="w")
        traj_writer = Trajectory(traj_file, "w", atoms)
        dyn.attach(logger, interval=dump)
        dyn.attach(traj_writer.write, interval=dump)
        dyn.run(nprod)

        frames = []
        for at in read(traj_file, ":"):
            forces = at.get_forces()
            displ = at.positions - atoms_ideal.get_positions()
            at.positions = atoms_ideal.get_positions()
            at.new_array("displacements", displ)
            at.new_array("forces", forces)
            frames.append(at.copy())
        print(f" Number of snapshots: {len(frames)}")
        out_xyz = f"md_runs/snapshots_T{T}.xyz"
        write(out_xyz, frames, format="extxyz")

        structures = read(out_xyz, index=":")
        cs = ClusterSpace(structures[0], [cutoffs[0]])
        stc = StructureContainer(cs)
        for s in structures:
            stc.add_structure(s)

        opt = Optimizer(stc.get_fit_data(), train_size=1.0)
        opt.train()
        print(opt)
        fcp_ehm = ForceConstantPotential(cs, opt.parameters)
        fcp_ehm.write(f"Potential_ehm_T{T}.fcp")

        fcs_ehm = fcp_ehm.get_force_constants(sc)
        fcs_ehm.write_to_phonopy(f"FORCE_CONSTANTS_2ND_{T}", format="text")


# -----------------------------------------------------------------------------#
# CLI
# -----------------------------------------------------------------------------#

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser()
    p.add_argument("-p", "--pc", type=str, default="pc.xyz",
                   help="Primitive cell structure (e.g., pc.xyz)")
    p.add_argument("-s", "--sc", type=str, default="sc.xyz",
                   help="Supercell structure (e.g., sc.xyz)")
    p.add_argument("-m", "--mode", type=str, default="pre",
                   help="Mode: pre, post, scph, ehm, shell, control")
    p.add_argument("-i", "--infiles", type=str, default="INFILES_CP2K",
                   help="Path to CP2K input template directory")
    p.add_argument("--coord_filename", type=str, default="coord.xyz",
                   help="Filename CP2K reads coordinates from")
    p.add_argument("--inp_glob", type=str, default="*.inp",
                   help="Glob to find CP2K input files to patch")
    p.add_argument("-r", "--rat_mode", type=str, default="rattle",
                   help="Distortion generator: rattle, mc, phon-rat")
    p.add_argument("-n", "--rat_n", type=int, default=20,
                   help="Number of structures to generate/read")
    p.add_argument("-a", "--rat_amp", type=float, default=0.02,
                   help="Rattle amplitude (std. dev., Å)")
    p.add_argument("-c", "--cutoffs", nargs="+", default=[5],
                   help="Cutoffs in Å or tokens: MAX, MAX2, MAX3 or "
                        "negative N for Nth neighbor shell")
    p.add_argument("--temp", nargs="+", default=[300],
                   help="Temperatures in K (SCPH/EHM/control)")
    p.add_argument("--scaleb", type=float, default=0.1,
                   help="scalebroad for CONTROL file")
    p.add_argument("--ngrid", nargs="+", default=["1", "1", "1"],
                   help="q-mesh density for CONTROL")
    p.add_argument("--scell", nargs="+", default=["1", "1", "1"],
                   help="Supercell size for CONTROL")
    p.add_argument("-f", "--fit_methods", nargs="+", default=["rfe"],
                   help="Fit method(s): least-squares, rfe, ardr, lasso, …")
    p.add_argument("-t", "--train_size", type=float, default=0.8,
                   help="Train fraction for model fitting")
    p.add_argument("-b", "--benchmark", action="store_true",
                   help="Sweep multiple train sizes")
    p.add_argument("--rot_sumrule", action="store_true",
                   help="Enforce rotational sum rules")
    p.add_argument("--scph_alpha", type=float, default=0.2,
                   help="SCPH mixing parameter alpha")
    p.add_argument("--scph_niter", type=int, default=50,
                   help="SCPH iterations")
    p.add_argument("--scph_nstr", type=int, default=100,
                   help="SCPH structures per iteration")
    p.add_argument("--ehm_neq_st", type=int, default=3000,
                   help="EHM MD equilibration steps")
    p.add_argument("--ehm_pro_st", type=int, default=3000,
                   help="EHM MD production steps")
    p.add_argument("--ehm_tim_st", type=float, default=1.0,
                   help="EHM MD time step (fs)")
    p.add_argument("--ehm_dump", type=int, default=100,
                   help="EHM snapshot/log interval")
    p.add_argument("--ehm_size", type=int, default=6,
                   help="EHM MD supercell size (N×N×N)")

    # New features you requested
    p.add_argument("--trim_iqr_k", type=float, default=0.0,
                   help="Outlier trimming: drop structures whose max |F| "
                        "exceeds Q3 + k·IQR (0 disables).")
    p.add_argument("--cutoff_sweep", nargs="+", default=None,
                   help="Optional sweep of the FIRST cutoff (Å), "
                        "e.g. --cutoff_sweep 4 6 8")
    p.add_argument("--bootstrap", type=int, default=0,
                   help="Extra shuffle-split validations per model.")
    p.add_argument("--report_json", type=str, default="fit_report.json",
                   help="Aggregate JSON report path.")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    pc = read(args.pc)
    sc = read(args.sc)

    if args.mode == "control":
        temp = float(args.temp[0])
        run_control(args.ngrid, args.scell, temp)
        return

    if args.mode == "shell":
        run_shell(sc)
        return

    if args.mode == "pre":
        run_pre(sc, args.rat_mode, args.rat_n, args.rat_amp,
                args.infiles, args.coord_filename, args.inp_glob)
        return

    cutoffs, folding, lim = build_cutoffs(sc, args.cutoffs)

    if args.mode == "post":
        run_post(
            pc=pc,
            sc=sc,
            cutoffs=cutoffs,
            folding=folding,
            lim=lim,
            rat_n=args.rat_n,
            fit_methods=[str(m) for m in args.fit_methods],
            train_fraction=float(args.train_size),
            benchmark=bool(args.benchmark),
            rot_sumrule=bool(args.rot_sumrule),
            trim_iqr_k=float(args.trim_iqr_k),
            cutoff_sweep=[float(x) for x in args.cutoff_sweep]
                         if args.cutoff_sweep else None,
            bootstrap=int(args.bootstrap),
            report_json=str(args.report_json),
            coord_filename=args.coord_filename,
        )
        return

    if args.mode == "scph":
        run_scph(
            sc=sc,
            cutoffs=cutoffs,
            temps=[float(t) for t in args.temp],
            alpha=args.scph_alpha,
            n_iterations=args.scph_niter,
            n_structures=args.scph_nstr,
            rot_sumrule=bool(args.rot_sumrule),
        )
        return

    if args.mode == "ehm":
        run_ehm(
            sc=sc,
            temps=[float(t) for t in args.temp],
            size=args.ehm_size,
            neq=args.ehm_neq_st,
            nprod=args.ehm_pro_st,
            dt_fs=args.ehm_tim_st,
            dump=args.ehm_dump,
            cutoffs=cutoffs,
        )
        return

    raise ValueError(f"Unknown mode: {args.mode!r}")


if __name__ == "__main__":
    main()
