#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import spglib
import numpy as np
import argparse
from ase import Atoms
from ase.io import read, write
from ase.build import make_supercell
from phonon_tools import get_xyz_files, get_restart_files  # reusing existing functions

# I/O Helper Functions

def read_xyz(file):
    return read(file, format="extxyz")

def read_cp2k(file):
    return read(file, format="cp2k-restart")

def write_xyz(file, atoms_obj):
    write(file, atoms_obj, format="extxyz")

# CP2K writing is unsupported; fall back to extxyz
def write_cp2k(file, atoms_obj):
    write(file, atoms_obj, format="extxyz")

# Main Functions

def create_symmetry_structures(file_format):
    # Choose file retrieval, reader, and writer functions based on format.
    if file_format == "xyz":
        files = get_xyz_files('.')
        read_func = read_xyz
        write_func = write_xyz
    else:
        files = get_restart_files('.')
        read_func = read_cp2k
        write_func = write_cp2k

    if not files:
        print(f"No {file_format} files found in the directory.")
        return

    structure = read_func(files[0])
    cell = (structure.get_cell(), structure.get_scaled_positions(), structure.get_atomic_numbers())
    precisions = [1e-2, 1e-3, 0.0004, 1e-4, 1e-5, 1e-6]

    with open('detected_space_groups.txt', 'w') as log_file:
        log_file.write("Symmetry Precision\tDetected Space Group\n")
        for symprec in precisions:
            std_cell = spglib.standardize_cell(cell, to_primitive=False, no_idealize=False, symprec=symprec)
            mod_structure = Atoms(
                symbols=std_cell[2],
                positions=np.dot(std_cell[1], std_cell[0]),
                cell=std_cell[0],
                pbc=True
            )
            sg = spglib.get_symmetry_dataset((std_cell[0], std_cell[1], std_cell[2]), symprec=0.00001)['international']
            log_file.write(f"{symprec}\t{sg}\n")
            filename = f'Zn_MOF5_H_{symprec}.xyz'
            write_func(filename, mod_structure)
            print(f"Modified structure saved as {filename}, detected space group: {sg}")

def check_symmetry(file_format, files=None):
    if file_format == "xyz":
        files = files or get_xyz_files('.')
        read_func = read_xyz
    else:
        files = files or get_restart_files('.')
        read_func = read_cp2k

    tolerances = [1e-2, 0.005, 1e-3, 1e-4, 0.0004, 1e-5, 1e-6]
    for file in files:
        try:
            atoms = read_func(file)
            for tol in tolerances:
                sg = spglib.get_symmetry_dataset(
                    (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers()),
                    symprec=tol
                )['international']
                print(f"{os.path.basename(file)} with tol {tol}: {sg}")
        except Exception as e:
            print(f"Error processing {file}: {e}")

def create_supercells(file_format, supercell_size=(2, 2, 2)):
    if file_format == "xyz":
        files = get_xyz_files('.')
        read_func = read_xyz
        write_func = write_xyz
    else:
        files = get_restart_files('.')
        read_func = read_cp2k
        write_func = write_cp2k

    # Check symmetry for all input files.
    check_symmetry(file_format, files)
    for file in files:
        atoms = read_func(file)
        supercell = make_supercell(atoms, [
            [supercell_size[0], 0, 0],
            [0, supercell_size[1], 0],
            [0, 0, supercell_size[2]]
        ])
        filename = f"SC_{supercell_size[0]}{supercell_size[1]}{supercell_size[2]}_{os.path.basename(file)}"
        write_func(filename, supercell)
        print(f"Supercell saved as {filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate symmetry-standardized structures, supercells, or check symmetry for input structures."
    )
    parser.add_argument("-m", "--mode", choices=["symmetry", "supercell", "check"], required=True,
                        help="Mode: symmetry, supercell, or check")
    parser.add_argument("-s", "--supercell", nargs=3, type=int, metavar=('X', 'Y', 'Z'),
                        default=(2, 2, 2),
                        help="Supercell size (e.g., -s 2 2 2)")
    parser.add_argument("-f", "--format", choices=["xyz", "cp2k"],
                        default="xyz",
                        help="Input/Output format: 'xyz' or 'cp2k'")
    args = parser.parse_args()

    if args.mode == "symmetry":
        create_symmetry_structures(args.format)
    elif args.mode == "supercell":
        if args.supercell:
            create_supercells(args.format, supercell_size=tuple(args.supercell))
        else:
            print("Error: Supercell size must be provided with -s when using supercell mode.")
    elif args.mode == "check":
        check_symmetry(args.format)
