#!/usr/bin/env python3

import argparse
import os
from aim2dat.io.phonopy import read_band_structure, read_total_density_of_states, read_atom_proj_density_of_states
from aim2dat.plots.band_structure_dos import BandStructureDOSPlot
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as plt

def get_phonopy_file(base_path):
    """
    Get the appropriate phonopy file, trying 'phonopy.yaml' first and falling back to 'phonopy_disp.yaml'.
    """
    phonopy_yaml = os.path.join(base_path, "phonopy.yaml")
    phonopy_disp_yaml = os.path.join(base_path, "phonopy_disp.yaml")

    if os.path.exists(phonopy_yaml):
        return phonopy_yaml
    elif os.path.exists(phonopy_disp_yaml):
        return phonopy_disp_yaml
    else:
        raise FileNotFoundError(f"Neither 'phonopy.yaml' nor 'phonopy_disp.yaml' found in {base_path}")


def load_band_structures(base_path):
    """
    Load band structures and reference cells for a given base path.

    Parameters:
        base_path (str): The base directory path containing the required input files.

    Returns:
        tuple: A tuple containing two lists:
            - band_structures (list): List of band structures.
            - ref_cells (list): List of reference cells.
    """
    # Initialize containers for results
    band_structures = []
    ref_cells = []

    # Define q-points and labels for the band structure path
    qpoints = [
        [
            [0, 0, 0],       # Γ
            [0.5, 0.0, 0.5], # X
            [0.5, 0.25, 0.75], # W
            [0.375, 0.375, 0.75], # K
            [0.0, 0.0, 0.0], # Γ
            [0.5, 0.5, 0.5], # L
            [0.625, 0.25, 0.625], # U
            [0.5, 0.25, 0.75], # W
            [0.5, 0.5, 0.5], # L
            [0.375, 0.375, 0.75], # K
        ]
    ]

    labels = ["Γ", "X", "W", "K", "Γ", "L", "U", "W", "L", "K"]

    # Retrieve the necessary file paths
    phonopy_file = get_phonopy_file(base_path)
    force_sets_file = os.path.join(base_path, "FORCE_SETS")

    # Read band structure and reference cell
    bs, rc = read_band_structure(
        phonopy_file,
        qpoints,
        npoints=255,
        force_sets_file_name=force_sets_file,
        path_labels=labels,
		phonopy_kwargs={"symprec": 0.1},
    )
    band_structures.append(bs)
    ref_cells.append(rc)

    return band_structures, ref_cells



def load_dos(base_path):
    """
    Load both TDOS and atom-projected DOS for a given base path.
    """
    pdos = []
    tdos = []

    phonopy_file = get_phonopy_file(base_path)
    file_path_forces = os.path.join(base_path, "FORCE_SETS")

    # Loading atom-projected DOS
    pdos.append(read_atom_proj_density_of_states(
        phonopy_file,
        force_sets_file_name=file_path_forces,
        mesh=20
    ))
    # Loading TDOS
    tdos.append(read_total_density_of_states(
        phonopy_file,
        force_sets_file_name=file_path_forces,
        mesh=20
    ))
    return pdos, tdos


def main():
    parser = argparse.ArgumentParser(description="Generate band structure and DOS plots.")
    parser.add_argument("--base-paths", nargs="+", required=True,
                        help="Base paths for methods (format: label:path).")
    parser.add_argument("--pdos-method", required=True,
                        help="The method label to use for projected DOS.")
    parser.add_argument("--output", default="combined_band_dos_plot.png",
                        help="Output file name for the plot.")
    args = parser.parse_args()

    # Process input arguments
    base_paths = {entry.split(":")[0]: entry.split(":")[1] for entry in args.base_paths}
    methods = list(base_paths.keys())
    pdos_method = args.pdos_method

    colors = ['#2980B9', '#66BB6A', 'orange', 'C1', 'C2']

    # Initialize the plot
    plot = BandStructureDOSPlot()
    plot.auto_set_axis_properties(set_y_label=False)
    plot.ratio = (10, 6)
    plot.subplot_ncols = 5
    plot.subplot_gridspec = [(0, 1, 0, 3), (0, 1, 3, 5)]
    plot.custom_linestyles = ["solid"]
    plot.custom_colors = colors
    plot.show_legend = True
    plot.x_label = [None, "HDOS (THz$^{-1}$)"]
    plot.y_label = ["Frequency (THz)", None]
    plot.subplot_adjust = {"top": 1, "right": 1}
    plot.y_range = [-1.0, 10.0]

    # Import band structures for all methods
    for method in methods:
        band_structures, ref_cells = load_band_structures(base_paths[method])
        plot.set_reference_cell(ref_cells[0])
        plot.import_band_structure(data_label=method, **band_structures[0])

    # Import projected DOS for the specified method
    pdos, tdos = load_dos(base_paths[pdos_method])
    plot.import_total_dos(data_label=pdos_method, **tdos[0])
    plot.import_projected_dos(
        pdos_method,
        pdos[0]["energy"],
        pdos[0]["pdos"],
        sum_kinds=True,
        sum_principal_qn=True,
        sum_magnetic_qn=True,
    )

    # Generate and save the plot
    final_plot = plot.plot(methods)
    final_plot.savefig(args.output, dpi=300)
    plt.close(final_plot)  # Ensure the plot is closed and not displayed
    print(f"Plot saved as {args.output}")


if __name__ == "__main__":
    main()
