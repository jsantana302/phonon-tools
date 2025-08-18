#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
from ase.io import read, write
from ase.build import make_supercell
from ase.visualize import view
import spglib

# Load structure (CIF, VASP POSCAR, etc.)
atoms = read('PBE.xyz')

# Apply symmetry detection
dataset = spglib.get_symmetry_dataset((atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers()))
print(f"Detected space group: {dataset['international']}")

# Create a 2x2x2 supercell
supercell = make_supercell(atoms, [[2,0,0],[0,2,0],[0,0,2]])

# Visualize
view(supercell)
