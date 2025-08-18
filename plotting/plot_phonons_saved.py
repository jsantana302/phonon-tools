#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
from aim2dat.io.phonopy import read_band_structure
from aim2dat.plots.band_structure_dos import BandStructurePlot

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

def get_force_file(base_path):
    """
    Get the appropriate force file, trying 'FORCE_SETS' first and falling back to 'FORCE_CONSTANTS'.
    """
    force_sets = os.path.join(base_path, "FORCE_SETS")
    force_constants = os.path.join(base_path, "FORCE_CONSTANTS")

    if os.path.exists(force_sets):
        return {"force_sets_file_name": "FORCE_SETS"}
    elif os.path.exists(force_constants):
        return {"force_constants_file_name": "FORCE_CONSTANTS"}
    else:
        raise FileNotFoundError(f"Neither 'FORCE_SETS' nor 'FORCE_CONSTANTS' found in {base_path}")

# Get the appropriate force file
try:
    force_file_name = get_force_file(os.getcwd())
except FileNotFoundError as e:
    raise FileNotFoundError(str(e))

# Define the band path and labels according to your specifications
band_path = [
    [0.0, 0.0, 0.0],  # Gamma
    [0.5, 0, 0.5],    # X
    [0.5, 0.25, 0.75], # W
    [0.375, 0.375, 0.75], # K
    [0.0, 0.0, 0.0],  # Gamma
    [0.5, 0.5, 0.5],  # L
    [0.625, 0.25, 0.625], # U
    [0.5, 0.25, 0.75], # W
    [0.5, 0.5, 0.5],  # L
    [0.375, 0.375, 0.75]     # K
]
path_labels = ["Γ", "X", "W", "K", "Γ", "L", "U", "W", "L", "K"]

# Use the get_phonopy_file function to determine the input file
try:
    input_file = get_phonopy_file(os.getcwd())
except FileNotFoundError as e:
    raise FileNotFoundError(str(e))

# Load the band structure data
band_structure, ref_cell = read_band_structure(
    input_file,
    [band_path],
    255,
    path_labels=path_labels,
    phonopy_kwargs={"symprec": 0.1},
    **force_file_name  # Pass either force_sets_file_name or force_constants_file_name
)

# Get the current folder name and set it as the plot title and PDF name
current_folder_name = os.path.basename(os.getcwd())

# Initialize the band structure plot
bands_plot = BandStructurePlot()
bands_plot.y_label = "Frequency in THz"
bands_plot.show_plot = False  # Set to False for saving without displaying
bands_plot.import_band_structure(data_label="test_band_structure", **band_structure)
bands_plot.y_range = [-1.0, 10.0]

# Generate the plot with the folder name as title
plot = bands_plot.plot("test_band_structure")
#plot = bands_plot.plot("test_band_structure", plot_title=f"{current_folder_name}")


# Save the plot with the folder name in the PDF filename
pdf_filename = f"ph_{current_folder_name}.pdf"
plot.savefig(pdf_filename)
