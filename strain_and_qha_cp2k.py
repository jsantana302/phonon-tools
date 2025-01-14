#!/usr/bin/env python3

import os
import shutil
import subprocess
import argparse
import numpy as np
from aim2dat.strct import StructureCollection
from aim2dat.strct.strct_manipulation import scale_unit_cell

print("Script started...", flush=True)

def run_cmd(cmd, timeout=300):
    """Runs a shell command with real-time output and error handling."""
    try:
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(timeout=timeout)
        if stdout:
            print(stdout.strip())
        if process.returncode != 0:
            print(f"Error executing command: {cmd}\n{stderr.strip()}")
    except subprocess.TimeoutExpired:
        print(f"Command timed out: {cmd}")
        process.kill()

def get_directory_list():
    """Returns a list of directories starting with 'comp' or 'exp', sorted appropriately."""
    print("Scanning for directories...", flush=True)
    directories = [d for d in os.listdir('.') if os.path.isdir(d) and (d.startswith('comp') or d.startswith('exp'))]
    comp_dirs = sorted([d for d in directories if d.startswith('comp')], key=lambda x: float(x.split('_')[1]), reverse=True)
    exp_dirs = sorted([d for d in directories if d.startswith('exp')], key=lambda x: float(x.split('_')[1]))
    ordered_dirs = comp_dirs + exp_dirs
    print(f"Found directories: {ordered_dirs}", flush=True)
    return ordered_dirs

def read_lattice_from_xyz(file_path):
    """Reads the lattice vectors from an XYZ file."""
    with open(file_path, 'r') as f:
        for line in f:
            if "Lattice=" in line:
                lattice_data = line.split('Lattice=')[1].split('"')[1]
                lattice = [float(x) for x in lattice_data.split()]
                if len(lattice) != 9:
                    raise ValueError(f"Invalid lattice data in {file_path}: {lattice}")
                return lattice
    raise ValueError(f"No Lattice field found in {file_path}")

def calculate_volume(lattice):
    """Calculates the volume of a lattice from its vectors."""
    lattice = np.array(lattice) if isinstance(lattice, list) else lattice
    matrix = lattice.reshape(3, 3)
    volume = np.abs(np.linalg.det(matrix))
    return volume

def create_job_scripts():
    """Copies all job scripts starting with 'cp2k_qha_' from the source directory into the current directory."""
    source_dir = "/user/j.santanaandreo/u12658/scripts/jobs/"
    cmd = f"cp {source_dir}cp2k_qha_* ."
    run_cmd(cmd)

def qha_pre():
    """Pre-processing mode: Organizes files, modifies inputs, and creates job scripts for batch processing."""
    files = os.listdir('.')
    comp_folders = []
    exp_folders = []
    original_dir = os.getcwd()  # Save the original working directory

    for file in files:
        if file.endswith('.xyz'):
            prefix = 'comp' if 'comp' in file else 'exp' if 'exp' in file else None
            if not prefix:
                continue
            value = file.split('_')[-1].replace('.xyz', '')
            folder_name = f"{prefix}_{value}"
            if not os.path.exists(folder_name):
                os.mkdir(folder_name)
            shutil.move(file, os.path.join(folder_name, file))

            # Find the first .inp file in the folder
            inp_files = [f for f in files if f.endswith('.inp')]
            if not inp_files:
                print("Error: No .inp files found in the current directory.")
                continue

            source_inp = os.path.join('.', inp_files[0])  # Take the first .inp file
            dest_inp = os.path.join(folder_name, inp_files[0])
            shutil.copy(source_inp, dest_inp)

            # Modify the .inp file
            with open(dest_inp, 'r') as inp_file:
                lines = inp_file.readlines()

            topology_start = None
            topology_end = None
            for i, line in enumerate(lines):
                if '&TOPOLOGY' in line:
                    topology_start = i
                elif topology_start is not None and line.strip().startswith('&END'):
                    topology_end = i
                    break

            if topology_start is None or topology_end is None:
                print(f"Error: TOPOLOGY section not found in {dest_inp}.")
                continue

            new_topology = (
                f"     &TOPOLOGY\n"
                f"         COORD_FILE_NAME {file}\n"
                f"         COORD_FILE_FORMAT XYZ\n"
                f"     &END TOPOLOGY\n"
            )

            lines = lines[:topology_start] + [new_topology] + lines[topology_end + 1:]

            with open(dest_inp, 'w') as inp_file:
                inp_file.writelines(lines)

            # Change into the folder and execute the script
            os.chdir(folder_name)
            try:
                run_cmd(f'~/scripts/cell_script.bash')
            finally:
                os.chdir(original_dir)  # Return to the original directory

            if prefix == 'comp':
                comp_folders.append(folder_name)
            elif prefix == 'exp':
                exp_folders.append(folder_name)

    create_job_scripts()
    print("Pre-processing completed successfully.")




def create_ev(e_v_dat_name, strain_values):
    """
    Creates two E_V.dat files:
    1. `e-v.dat`: Contains only volume and energy columns.
    2. `e-v-full.dat`: Contains volume, energy, and directory name.
    """
    directories = get_directory_list()
    e_v_dat_full_name = e_v_dat_name.replace('.dat', '-full.dat')

    with open(e_v_dat_name, 'w') as ev_file, open(e_v_dat_full_name, 'w') as ev_full_file:
        ev_file.write("# Cell volume [Å³]   Total energy [a.u.]\n")
        ev_full_file.write("# Cell volume [Å³]   Total energy [a.u.]   Directory name\n")

        for dir in directories:
            out_files = [f for f in os.listdir(dir) if f.endswith('.out') and not f.startswith('slurm')]
            if len(out_files) == 0:
                print(f"Warning: No CP2K .out file found in {dir}.")
                continue

            out_file_path = os.path.join(dir, out_files[0])
            volume, energy = None, None

            with open(out_file_path, 'r') as f:
                for line in f:
                    if "CELL| Volume [angstrom^3]:" in line:
                        volume = float(line.split()[-1])
                    if "ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:" in line:
                        energy = float(line.split()[-1])
                    if volume is not None and energy is not None:
                        break

            if volume is not None and energy is not None:
                ev_file.write(f"{volume:<25} {energy:<25}\n")
                ev_full_file.write(f"{volume:<25} {energy:<25} {dir}\n")
            else:
                print(f"Warning: Volume or energy not found in {out_file_path}.")

    print(f"Files created: {e_v_dat_name} and {e_v_dat_full_name}.")


def apply_strain(volume_increment, num_volumes):
    """Applies volume strain with a uniform volume increment specified as arguments."""
    xyz_files = [f for f in os.listdir('.') if f.endswith('.xyz')]
    if not xyz_files:
        print("Error: No XYZ files found in the current directory.")
        return

    input_file = xyz_files[0]
    print(f"Using input file: {input_file}")

    # Extract lattice and volume
    lattice = read_lattice_from_xyz(input_file)
    original_volume = calculate_volume(lattice)

    # Generate target volumes
    target_volumes = [
        original_volume + i * volume_increment for i in range(-num_volumes, num_volumes + 1)
    ]

    # Load structure
    structure_name = os.path.splitext(os.path.basename(input_file))[0]
    strct = StructureCollection()
    strct.append_from_file(structure_name, input_file)
    original_structure = strct[0]

    # Apply strains
    for i, target_volume in enumerate(target_volumes):
        scaling_factor = (target_volume / original_volume) ** (1 / 3)
        strain_type = "exp" if target_volume > original_volume else "comp"
        strain_value = i - num_volumes

        # Generate new structure with scaled volume
        new_structure = scale_unit_cell(original_structure, scaling_factors=scaling_factor, change_label=True)
        new_filename = f"{structure_name}_{strain_type}_{abs(int(strain_value))}.xyz"
        new_structure.to_file(new_filename)

    print("Strain application completed successfully.")

def submit_phonopy():
    """
    Submits phonopy calculations for each subdirectory using the specified submit command.
    """
    print("Submitting phonopy calculations...")
    main_directories = get_directory_list()
    if not main_directories:
        print("No main directories found for submission.")
        return

    submit_command = os.path.expanduser("~/scripts/nodes_and_send.sh")

    for main_dir in main_directories:
        try:
            # Find the subdirectory with the naming convention
            subdirectory_name = f"pbe_D3_222_prim_{os.path.basename(main_dir)}"
            subdirectory_path = os.path.join(main_dir, subdirectory_name)

            if not os.path.exists(subdirectory_path):
                print(f"Subdirectory {subdirectory_path} does not exist. Skipping.")
                continue

            # Run the submission command inside the subdirectory
            print(f"Submitting job for {subdirectory_path} using {submit_command}")
            run_cmd(f"cd {subdirectory_path} && {submit_command}")

        except Exception as e:
            print(f"Error processing {main_dir}: {e}")

def phonopy_pre(tolerance):
    """
    Prepares subdirectories for phonopy post-processing by ensuring folders,
    copying required files, and running the preparation script.
    """
    print("Starting phonopy pre-processing...")
    main_directories = get_directory_list()
    if not main_directories:
        print("No main directories found for pre-processing.")
        return

    for main_dir in main_directories:
        # Find the restart file in the current main directory
        restart_files = [f for f in os.listdir(main_dir) if f.endswith('-1.restart')]
        if not restart_files:
            print(f"No restart file found in {main_dir}. Skipping.")
            continue
        restart_file = restart_files[0]
        base_name = restart_file.split('-1.restart')[0]

        # Create or use the subdirectory with the desired naming convention
        subdirectory_name = f"pbe_D3_222_prim_{os.path.basename(main_dir)}"
        subdirectory_path = os.path.join(main_dir, subdirectory_name)
        print(f"Processing subdirectory: {subdirectory_path}")

        try:
            if not os.path.exists(subdirectory_path):
                os.makedirs(subdirectory_path)
                print(f"Created subdirectory: {subdirectory_path}")
            else:
                print(f"Subdirectory already exists: {subdirectory_path}")

            # Copy the restart file and rename it if not already present
            destination_input_file = os.path.join(subdirectory_path, f"{base_name}.inp")
            if not os.path.exists(destination_input_file):
                shutil.copy(os.path.join(main_dir, restart_file), destination_input_file)
                print(f"Copied {restart_file} to {destination_input_file}")
            else:
                print(f"Input file already exists: {destination_input_file}")

            # Execute preparation script inside the subdirectory
            preparation_script = os.path.expanduser("~/scripts/pcp2k.py")
            preparation_command = f"python {preparation_script} -m pre -t {tolerance}"
            run_cmd(f"cd {subdirectory_path} && {preparation_command}")

        except Exception as e:
            print(f"Error processing {main_dir}: {e}")

def phonopy_post():
    """Executes the phonopy post-processing script inside each relevant folder."""
    print("Starting phonopy post-processing...")
    directories = get_directory_list()
    if not directories:
        print("No directories found for processing.")
        return

    initial_dir = os.getcwd()
    for dir in directories:
        subfolders = [os.path.join(dir, f) for f in os.listdir(dir) if f.startswith('pbe_D3_222_prim_') and os.path.isdir(os.path.join(dir, f))]
        if not subfolders:
            print(f"No subfolders found in {dir}. Skipping.")
            continue

        for subfolder in subfolders:
            print(f"Processing: {subfolder}")
            try:
                os.chdir(subfolder)
                script_path = os.path.expanduser("~/scripts/pcp2k.py")
                log_file = f"post_{os.path.basename(subfolder)}.log"
                command = f"python {script_path} -m post -tp > {log_file} 2>&1 &"
                print(f"Executing command: {command}")
                os.system(command)
                print(f"Post-processing started in the background for {subfolder}")
            except Exception as e:
                print(f"Error processing {subfolder}: {e}")
            finally:
                os.chdir(initial_dir)


def qha_post(temp, pressure):
    """
    Performs QHA analysis by finding and processing all thermal_properties.yaml files
    from relevant directories provided by get_directory_list. Copies all found files to
    a new folder named 'thermal_properties' with renamed files based on their parent folder.
    """

    print("Starting QHA analysis...")

    # Use get_directory_list to get the ordered directories
    directories = get_directory_list()

    # Ensure the 'thermal_properties' folder exists
    thermal_properties_dir = "thermal_properties"
    os.makedirs(thermal_properties_dir, exist_ok=True)

    # Search for thermal_properties.yaml files in the subdirectories of the listed directories
    thermal_property_files = []
    for directory in directories:
        for root, _, files in os.walk(directory):
            for file in files:
                if file == "thermal_properties.yaml":
                    file_path = os.path.join(root, file)
                    thermal_property_files.append(file_path)
                    # Rename file with folder name and copy to the thermal_properties folder
                    parent_folder = os.path.basename(root)
                    new_file_name = f"thermal_properties_{parent_folder}.yaml"
                    shutil.copy(file_path, os.path.join(thermal_properties_dir, new_file_name))

    if not thermal_property_files:
        print("No thermal_properties.yaml files found. Exiting QHA analysis.")
        return

    # Generate the phonopy-qha command with the located files
    cmd = f"phonopy-qha -p -s --tmax={temp} e-v.dat " + " ".join(thermal_property_files)
    print(f"Executing: {cmd}")
    run_cmd(cmd)


def main():
    parser = argparse.ArgumentParser(description="Script with multiple modes for processing structures.")
    parser.add_argument(
        '-m', '--mode', required=True,
        choices=['qha_pre', 'create_ev', 'apply_strain', 'phonopy_pre', 'phonopy_post', 'qha_post', 'submit_phonopy'],
        help="Mode of operation."
    )
    parser.add_argument('-v', '--volume', type=float, default=20,
                        help='Volume increment for strain (default: 20 Å³).')
    parser.add_argument('-nv', '--number_volumes', type=int, default=10,
                        help='Number of volumes to generate on each side of original (default: 4).')
    parser.add_argument('--e_vdat_name', default="e-v.dat", help="Name of the output E_V.dat file (default: e-v.dat).")
    parser.add_argument('-t', '--tolerance', type=float, default=0.04,
                        help='Tolerance value for phonopy_pre (default: 0.04).')
    parser.add_argument('--temp', default=1000, type=int, help='Temperature for QHA (default: 1000 K).')
    parser.add_argument('--pressure', default=0, type=int, help='Pressure for QHA (default: 0 GPa).')

    args = parser.parse_args()

    if args.mode == 'qha_pre':
        qha_pre()
    elif args.mode == 'create_ev':
        create_ev(args.e_vdat_name, [])
    elif args.mode == 'apply_strain':
        apply_strain(volume_increment=args.volume, num_volumes=args.number_volumes)
    elif args.mode == 'phonopy_pre':
        phonopy_pre(args.tolerance)
    elif args.mode == 'phonopy_post':
        phonopy_post()
    elif args.mode == 'qha_post':
        qha_post(args.temp, args.pressure)
    elif args.mode == 'submit_phonopy':
        submit_phonopy()
    else:
        print(f"Error: Unknown mode {args.mode}.")
        sys.exit(1)

if __name__ == "__main__":
    main()
