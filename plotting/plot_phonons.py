#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
from aim2dat.io.phonopy import read_band_structure
from aim2dat.plots.band_structure_dos import BandStructurePlot
from phonon_tools import get_phonopy_file, get_force_file  

force_file_name = get_force_file(os.getcwd())

band_path = [
    [0.0, 0.0, 0.0],
    [0.5, 0, 0.5],
    [0.5, 0.25, 0.75],
    [0.375, 0.375, 0.75],
    [0.0, 0.0, 0.0],
    [0.5, 0.5, 0.5],
    [0.625, 0.25, 0.625],
    [0.5, 0.25, 0.75],
    [0.5, 0.5, 0.5],
    [0.375, 0.375, 0.75]
]
path_labels = ["Γ", "X", "W", "K", "Γ", "L", "U", "W", "L", "K"]

input_file = get_phonopy_file(os.getcwd())

band_structure, ref_cell = read_band_structure(
    input_file,
    [band_path],
    255,
    path_labels=path_labels,
    phonopy_kwargs={"symprec": 0.1},
    **force_file_name
)

current_folder_name = os.path.basename(os.getcwd())

# Determine color based on folder name prefix

if current_folder_name.startswith("r2scan_"):
    plot_color = "C1" #orange
    plot_y_range= [-1.0, 8.0]
elif current_folder_name.startswith("pbe0_"):
    plot_color = "C3" #red
    plot_y_range = [-1.0, 8.0]
else:
    plot_color = "C0"  #blue
    plot_y_range = [-1.0, 8.0]


bands_plot = BandStructurePlot()
bands_plot.y_label = " "
#bands_plot.y_label = "Frequency in THz"
bands_plot.show_plot = False
bands_plot.import_band_structure(data_label="test_band_structure", **band_structure)
bands_plot.y_range = plot_y_range
bands_plot.custom_colors = [plot_color]  # Set custom colors based on folder name

plot = bands_plot.plot("test_band_structure")

pdf_filename = f"ph_{current_folder_name}.pdf"
jpg_filename = f"ph_{current_folder_name}.jpg"

plot.savefig(pdf_filename)

# Remove y-label for JPG
#bands_plot.y_label = " "
plot = bands_plot.plot("test_band_structure")
plot.savefig(jpg_filename, dpi=600, bbox_inches='tight')  # High DPI and no title
