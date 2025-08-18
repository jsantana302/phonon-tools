#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
from aim2dat.io.phonopy import read_total_density_of_states, read_atom_proj_density_of_states
from aim2dat.plots.band_structure_dos import DOSPlot

def filter_data(energy, dos, cutoff):
    """Filter the energy and DOS lists to energies <= cutoff."""
    filtered_energy = []
    filtered_dos = []
    for e, d in zip(energy, dos):
        if e <= cutoff:
            filtered_energy.append(e)
            filtered_dos.append(d)
    return filtered_energy, filtered_dos

def filter_projected_dos(energy, pdos, cutoff):
    """Filter the projected DOS data. `pdos` is expected to be a list of dictionaries."""
    mask = [e <= cutoff for e in energy]
    filtered_energy = [e for e in energy if e <= cutoff]
    filtered_pdos = []
    for dos_dict in pdos:
        new_dict = {}
        for key, values in dos_dict.items():
            new_values = [v for v, m in zip(values, mask) if m]
            new_dict[key] = new_values
        filtered_pdos.append(new_dict)
    return filtered_energy, filtered_pdos

# Read the DOS data
pdos = read_atom_proj_density_of_states(
    "phonopy.yaml",
    force_sets_file_name="FORCE_SETS",
    mesh=100,
)
tdos = read_total_density_of_states(
    "phonopy.yaml",
    force_sets_file_name="FORCE_SETS",
    mesh=100,
)

# Set the cutoff (20 THz)
cutoff = 20

# Filter the total DOS data
filtered_energy_total, filtered_tdos = filter_data(tdos["energy"], tdos["tdos"], cutoff)
# Reconstruct the total DOS dictionary (preserving unit_x if needed)
filtered_tdos_dict = {
    "energy": filtered_energy_total,
    "tdos": filtered_tdos,
    "unit_x": tdos.get("unit_x", "THz")
}

# Filter the projected DOS data
filtered_energy_proj, filtered_pdos = filter_projected_dos(pdos["energy"], pdos["pdos"], cutoff)

# Set up the DOS plot
dos_plot = DOSPlot()
dos_plot.y_label = "DOS in states/THz/cell"
dos_plot.x_label = "Frequency (THz)"
dos_plot.y_range = [-1, 100]
dos_plot.x_range = [0, 50]
dos_plot.show_legend = True

# Import filtered DOS data
dos_plot.import_projected_dos(
    "test_dos",
    filtered_energy_proj,
    filtered_pdos,
    sum_kinds=True,
    sum_principal_qn=True,
    sum_magnetic_qn=True,
)
dos_plot.import_total_dos("test_dos", **filtered_tdos_dict)

plot = dos_plot.plot("test_dos")
plot.savefig("Dos_filtered.png", dpi=300)

