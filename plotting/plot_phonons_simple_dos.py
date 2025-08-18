#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
from aim2dat.io.phonopy import (
    read_band_structure,
    read_total_density_of_states,
    read_atom_proj_density_of_states,
)
from aim2dat.plots.band_structure_dos import BandStructureDOSPlot

PHONOPY_KW = {
    "primitive_matrix": "auto",
    "symprec": 1e-3,
    }
def load_band_structure(current_dir="."):
    qpoints = [[
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.25, 0.75],
        [0.375, 0.375, 0.75],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.625, 0.25, 0.625],
        [0.5, 0.25, 0.75],
        [0.5, 0.5, 0.5],
        [0.375, 0.375, 0.75],
    ]]
    labels = ["Γ", "X", "W", "K", "Γ", "L", "U", "W", "L", "K"]
    phonopy_file = os.path.join(current_dir, "phonopy_disp.yaml")
    force_sets_file = os.path.join(current_dir, "FORCE_SETS")
    bs, rc = read_band_structure(
        phonopy_file,
        qpoints,
        npoints=51,
        force_sets_file_name=force_sets_file,
        path_labels=labels,
        phonopy_kwargs=PHONOPY_KW,
    )
    return bs, rc

def load_total_dos(current_dir="."):
    phonopy_file = os.path.join(current_dir, "phonopy_disp.yaml")
    force_sets_file = os.path.join(current_dir, "FORCE_SETS")
    tdos = read_total_density_of_states(
        phonopy_file,
        force_sets_file_name=force_sets_file,
        phonopy_kwargs=PHONOPY_KW,
    )
    return tdos

def load_projected_dos(current_dir="."):
    phonopy_file = os.path.join(current_dir, "phonopy_disp.yaml")
    force_sets_file = os.path.join(current_dir, "FORCE_SETS")
    pdos = read_atom_proj_density_of_states(
        phonopy_file,
        force_sets_file_name=force_sets_file,
        phonopy_kwargs=PHONOPY_KW,
    )
    return pdos

def main():
    base_path = "."
    band_structure, ref_cell = load_band_structure(base_path)
    total_dos = load_total_dos(base_path)
    pdos = load_projected_dos(base_path)

    colors = ["C0", "C1", "C2", "C3", "C4"]
    plot = BandStructureDOSPlot()
    plot.show_plot = False
    plot.store_plot = True
    plot.auto_set_axis_properties(set_y_label=False)
    plot.ratio = (10, 6)
    plot.subplot_ncols = 5
    plot.subplot_gridspec = [(0, 1, 0, 3), (0, 1, 3, 5)]
    plot.custom_linestyles = ["solid"]
    plot.custom_colors = colors
    plot.x_label = [None, r"DOS (THz$^{-1}$)"]
    plot.y_label = ["Frequency (THz)", None]
    plot.subplot_adjust = {"top": 1, "right": 1}
    plot.y_range = [-1.0, 4.5]
    plot.x_range = [0.0, 100.0]
    plot.show_legend = [False, True]

    plot.set_reference_cell(ref_cell)
    plot.import_band_structure(data_label="Phonon Dispersion",
                               **band_structure)
    plot.import_total_dos(data_label="Phonon Dispersion", **total_dos)
    plot.import_projected_dos(
        "Phonon Dispersion",
        pdos["energy"],
        pdos["pdos"],
        sum_kinds=True,
        sum_principal_qn=True,
        sum_magnetic_qn=True,
    )
    plot_name = "phonons_dos_highsymprec.png"
    _ = plot.plot(
        "Phonon Dispersion",
        plot_title="",
        plot_name=plot_name
    )
    print("Plot saved as " + plot_name)

if __name__ == "__main__":
    main()
