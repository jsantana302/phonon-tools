#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
from aim2dat.io.phonopy import read_total_density_of_states
from aim2dat.plots.band_structure_dos import DOSPlot
from aim2dat.io.phonopy import read_atom_proj_density_of_states

pdos = read_atom_proj_density_of_states(
    "phonopy.yaml",
    force_sets_file_name="FORCE_SETS",
    mesh=50,
)
tdos = read_total_density_of_states(
    "phonopy.yaml",
    force_sets_file_name="FORCE_SETS",
    mesh=50,
)

dos_plot = DOSPlot()
dos_plot.y_label = "DOS in states/THz/cell"
dos_plot.x_label = "Frequency (THz)"
dos_plot.y_range = [-1, 100]
dos_plot.x_range = [ 0.0, 4.5]
dos_plot.show_legend = True

dos_plot.import_projected_dos(
    "test_dos",
    pdos["energy"],
    pdos["pdos"],
    sum_kinds=True,
    sum_principal_qn=True,
    sum_magnetic_qn=True,
)


dos_plot.import_total_dos("test_dos", **tdos)
plot = dos_plot.plot("test_dos")
plot.savefig("Dos_frec_4.5.png", dpi=300)
