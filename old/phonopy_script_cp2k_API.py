#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import shutil
import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import write_supercells_with_displacements
from phonopy.cui import show_symmetry


# Function to run shell commands
def run_cmd(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error executing command: {cmd}\n{result.stderr}")
    return result.stdout

# Path to the Conda environment's Python interpreter
conda_py = "/home/nipj5103/miniconda3/envs/phonopyenv/bin/python"

# Variables
task = "pre"  # pre, post, or qha
mat = "Si"  # Material name

# Pre-processing variables
cell_type = "prim"  # prim or conv
sc = [[2, 0, 0], [0, 2, 0], [0, 0, 2]] # Supercell sizes


# QHA variables
n_files = 9
temp = 980
pressure = 0

# Directory for CP2K files
cp2k_dir = f"/home/nipj5103/files_CP2K/{mat}"

#Read crystal structure
unitcell, optional_structure_info = read_crystal_structure("{mat}.inp", interface_mode='cp2k')


# Perform operations for pre-processing
if task == "pre":
    # Symmetry and tolerance
    check_symmetry(unitcell, optional_structure_info)

    # Copy the appropriate unit cell file
    if cell_type == "prim":
        shutil.copy("Punitcell.inp", f"{mat}.inp")
    elif cell_type == "conv":
        shutil.copy("Bunitcell.inp", f"{mat}.inp")
    else:
        print("I don't understand")
        exit()

    # Create the supercell
	phonon = Phonopy(unitcell,
                 supercell_matrix= sc,
                 primitive_matrix=[[0, 0.5, 0.5],
                                   [0.5, 0, 0.5],
                                   [0.5, 0.5, 0]])	
    phonopy.generate_displacements(distance=0.01)
    supercells = phonon.supercells_with_displacements
	write_supercells_with_displacements("{mat}.in", supercell, interface_mode='cp2k', optional_structure_info=optional_structure_info)


# Perform operations for post-processing
elif task == "post":
    # Count the number of supercell files
    supercell_files = [f for f in os.listdir() if f.startswith(f"{mat}-supercell-") and f.endswith("-forces-1_0.xyz")]
    forces = [read_cp2k(os.path.join(os.getcwd(), filename)) for filename in supercell_files]
    
    # Read the force constants
    phonopy = Phonopy(get_primitive(f"{mat}.inp"), supercell_matrix=[s1, s2, s3])
    phonopy.produce_force_constants(forces)

    # Save force constants and phonon band structure
    phonopy.save()
    phonopy.write_yaml()

    # Create and save band structure
    phonopy.run_band_structure()
    phonopy.plot_band_structure().show()

# Perform operations for QHA
elif task == "qha":
    from phonopy.qha import PhonopyQHA

    # Prepare QHA input files
    thermal_properties_files = [f"thermal_properties-{i:03d}.yaml" for i in range(1, n_files + 1)]
    phonopy_qha = PhonopyQHA(thermal_properties_files, tmax=temp, pressure=pressure)

    # Save QHA results
    phonopy_qha.write_yaml()
    phonopy_qha.plot_qha().show()

else:
    print("I don't understand")
