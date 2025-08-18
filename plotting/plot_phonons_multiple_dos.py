#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from aim2dat.io.phonopy import (
    read_band_structure,
    read_total_density_of_states,
    read_atom_proj_density_of_states,
)
from aim2dat.plots.band_structure_dos import BandStructureDOSPlot
from phonon_tools import get_phonopy_file

def load_band_structures(bp):
    q = [[
        [0, 0, 0], [0.5, 0, 0.5], [0.5, 0.25, 0.75],
        [0.375, 0.375, 0.75], [0, 0, 0], [0.5, 0.5, 0.5],
        [0.625, 0.25, 0.625], [0.5, 0.25, 0.75],
        [0.5, 0.5, 0.5], [0.375, 0.375, 0.75],
    ]]
    lbl = ["Γ", "X", "W", "K", "Γ", "L", "U", "W", "L", "K"]
    pf = get_phonopy_file(bp)
    fs = os.path.join(bp, "FORCE_SETS")
    bs, rc = read_band_structure(
        pf, q,
        npoints=255,
        force_sets_file_name=fs,
        path_labels=lbl,
        phonopy_kwargs={"symprec": 1e-5},
    )
    return [bs], [rc]

def load_dos(bp):
    disp = os.path.join(bp, "phonopy_disp.yaml")
    fs = os.path.join(bp, "FORCE_SETS")
    if not os.path.isfile(disp):
        raise FileNotFoundError(f"{disp} not found")
    if not os.path.isfile(fs):
        raise FileNotFoundError(f"{fs} not found")
    tdos = read_total_density_of_states(disp, force_sets_file_name=fs)
    pdos = read_atom_proj_density_of_states(disp, force_sets_file_name=fs)
    return [pdos], [tdos]

def main():
    p = argparse.ArgumentParser(
        description="Multi phonon+PDOS plot"
    )
    p.add_argument(
        "--base-paths", nargs="+", required=True,
        help="label:path"
    )
    p.add_argument(
        "--pdos-method", required=True,
        help="which label to use for PDOS"
    )
    p.add_argument(
        "--output", default="combined.png",
        help="output filename"
    )
    args = p.parse_args()

    bp = dict(e.split(":", 1) for e in args.base_paths)
    methods = list(bp.keys())
    pdm = args.pdos_method

    label_map = {
        "eq": "eq",
        "comp_3": "1.3%",
        "comp_5": "2.2%",
        "comp_6": "2.6%",
    }
    display_labels = [label_map.get(m, m) for m in methods]

    plot = BandStructureDOSPlot()
    plot.auto_set_axis_properties(set_y_label=False)
    plot.ratio = (10, 6)
    plot.subplot_ncols = 5
    plot.subplot_gridspec = [(0, 1, 0, 3), (0, 1, 3, 5)]
    plot.custom_linestyles = ["solid"]
    plot.custom_colors = ["C0", "C3", "C2", "C1", "C4","C5", "C6", "orange"]
    plot.x_label = [None, "DOS (THz$^{-1}$)"]
    plot.y_label = ["Frequency (THz)", None]
    plot.subplot_adjust = {"top": 1, "right": 1}
    plot.y_range = [-1.0, 3.0]
    plot.x_range = [0, 100.0]
    plot.show_legend = True

    for m in methods:
        disp = label_map.get(m, m)
        bs, rc = load_band_structures(bp[m])
        plot.set_reference_cell(rc[0])
        plot.import_band_structure(data_label=disp, **bs[0])

    pdl, tdl = load_dos(bp[pdm])
    disp = label_map.get(pdm, pdm)
    plot.import_total_dos(data_label=disp, **tdl[0])
    plot.import_projected_dos(
        disp,
        pdl[0]["energy"],
        pdl[0]["pdos"],
        sum_kinds=True,
        sum_principal_qn=True,
        sum_magnetic_qn=True,
    )

    fig = plot.plot(display_labels)
    fig.savefig(args.output, dpi=300)
    plt.close(fig)

if __name__ == "__main__":
    main()
