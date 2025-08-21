#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import shutil
import argparse
import csv
import glob
import sys
import numpy as np
from aim2dat.strct import Structure, StructureCollection
from aim2dat.strct.strct_manipulation import scale_unit_cell
from pcp2k import pre_mode as pcp2k_pre
from phonon_tools import (
    get_directory_list,
    read_lattice_from_xyz,
    run_cmd,
    get_xyz_files,
    get_inp_files,
    get_cif_files,
    get_restart_files,
    erase_coord
)


print("Script started...", flush=True)


def calculate_volume(lattice):
    """Calculates the volume of a lattice from its vectors."""
    lattice = np.array(lattice) if isinstance(lattice, list) else lattice
    matrix = lattice.reshape(3, 3)
    volume = np.abs(np.linalg.det(matrix))
    return volume


def apply_strain(volume_increment, num_volumes):
        
    # Get the input structure file
    xyz_files = [f for f in os.listdir('.') if f.endswith('.xyz')]
    if not xyz_files:
        print("Error: No XYZ files found.")
        return
    input_file = xyz_files[0]
    print(f"Using input file: {input_file}")

    # Read structure and compute original volume
    structure_name = os.path.splitext(input_file)[0]
    lattice = read_lattice_from_xyz(input_file)
    original_volume = calculate_volume(lattice)

    # Create target volumes
    target_volumes = [original_volume + i * volume_increment for i in range(-num_volumes, num_volumes + 1)]
    
    # Load original structure
    strct = StructureCollection()
    strct.append_from_file(structure_name, input_file)
    original_structure = strct[0]

    # Apply strains
    for i, target_volume in enumerate(target_volumes):
        scaling_factor = (target_volume / original_volume) ** (1 / 3)
        strain_type = "exp" if target_volume > original_volume else "comp"
        strain_value = i - num_volumes

        # Scale the structure
        new_structure_dict = scale_unit_cell(original_structure, scaling_factors=scaling_factor, change_label=True)
        
        # Convert dictionary to Structure
        new_structure = Structure(**new_structure_dict)

        # Save to file
        new_filename = f"{structure_name}_{strain_type}_{abs(int(strain_value))}.xyz"
        new_structure.to_file(new_filename)

    print("Strain application completed successfully.")

def qha_pre():
    """
    Pre-processing mode: Organizes files, modifies inputs, and creates job
    scripts for batch processing.
    """
    xyz_files = get_xyz_files('.')  # Use the new function to get sorted .xyz files
    if not xyz_files:
        print("Error: No .xyz files found in the current directory.")
        return  # Exit early if no XYZ files are found

    comp_folders = []
    exp_folders = []
    original_dir = os.getcwd()
    job_script_path = "/user/j.santanaandreo/u12658/scripts/jobs"

    # Get all cp2k_geo* files in the job_script_path
    geo_files = glob.glob(os.path.join(job_script_path, 'cp2k_geo_*'))

    # Copy all cp2k_geo* files to the parent directory (current working directory)
    for geo_file in geo_files:
        shutil.copy(geo_file, os.path.basename(geo_file))

    for file in xyz_files:  # Iterate over the sorted .xyz files
        prefix = (
            'comp' if 'comp' in file else
            'exp' if 'exp' in file else
            None
        )
        if not prefix:
            continue

        value = file.split('_')[-1].replace('.xyz', '')
        folder_name = f"{prefix}_{value}"
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)

        shutil.move(file, os.path.join(folder_name, file))

        inp_files = get_inp_files('.')
        if not inp_files:
            continue  # Skip this iteration if no .inp files found

        source_inp = os.path.join('.', inp_files[0])
        dest_inp = os.path.join(folder_name, inp_files[0])
        shutil.copy(source_inp, dest_inp)

        # Modify the input file erasing the COORD section
        erase_coord(dest_inp)

        # Modify input file topology section
        with open(dest_inp, 'r') as inp_file:
            lines = inp_file.readlines()

        topology_start, topology_end = None, None
        for i, line in enumerate(lines):
            if "&TOPOLOGY" in line:
                topology_start = i
            elif topology_start is not None and line.strip().startswith("&END"):
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

        # Run command inside the newly created folder
        os.chdir(folder_name)
        try:
            run_cmd("~/scripts/cell_script.bash")
        finally:
            os.chdir(original_dir)

        # Store processed folders
        if prefix == 'comp':
            comp_folders.append(folder_name)
        elif prefix == 'exp':
            exp_folders.append(folder_name)

    print("Pre-processing completed successfully.")


def submit_geo():
    """
    Submits geometry optimization jobs for each directory using the specified
    submission command. Copies the job script to each directory before submission.
    """
    directories = get_directory_list(None)
    job_script_path = "/user/j.santanaandreo/u12658/scripts/jobs/cp2k_geo.job"
    
    if not os.path.exists(job_script_path):
        print(f"Job script {job_script_path} not found. Aborting.")
        return
    
    for directory in directories:
        job_script_dest = os.path.join(directory, os.path.basename(job_script_path))
        
        # Copy the job script to the directory
        shutil.copy(job_script_path, job_script_dest)
        print(f"Copied {job_script_path} to {directory}.")
        
        # Submit the job
        print(f"Submitting job for {directory} using {os.path.basename(job_script_path)}")
        run_cmd(f"cd {directory} && sbatch {os.path.basename(job_script_path)}")

def create_ev():
    e_v_dat = "e-v.dat"
    e_v_csv = "e-v-full.csv"
    dirs = get_directory_list()
    au_to_ev = 27.211386245988

    with open(e_v_dat, "w") as datf, open(e_v_csv, "w", newline="") as csvf:
        writer = csv.writer(csvf)
        datf.write("# Cell volume [Å³]   Total energy [eV]\n")
        writer.writerow(["Cell volume [Å³]", "Total energy [eV]", "Directory"])

        for d in dirs:
            outs = [f for f in os.listdir(d) if f.endswith(".out")
                    and not f.startswith("slurm")]
            if not outs:
                print(f"Warning: no .out in {d}")
                continue

            path = os.path.join(d, outs[0])
            last_vol = last_en = None
            with open(path) as f:
                for line in f:
                    if "CELL| Volume [angstrom^3]:" in line:
                        last_vol = float(line.split()[-1])
                    if ("ENERGY| Total FORCE_EVAL ( QS ) energy"
                            " [a.u.]:") in line:
                        last_en = float(line.split()[-1]) * au_to_ev

            if last_vol is not None and last_en is not None:
                datf.write(f"{last_vol:<25}{last_en:<25}\n")
                writer.writerow([last_vol, last_en, d])
            else:
                print(f"Warning: incomplete data in {path}")

    print(f"Created {e_v_dat} and {e_v_csv}")


def phonopy_pre(tolerance, cell_tipe, supercell_size):

    print("Starting phonopy pre-processing...")
    main_directories = get_directory_list()
    if not main_directories:
        print("No main directories found for pre-processing.")
        return

    supercell_str = "".join(map(str, supercell_size))

    for main_dir in main_directories:
        restart_files = get_restart_files(main_dir)
        if not restart_files:
            print(f"No restart file found in {main_dir}. Skipping.")
            continue
        restart_file = restart_files[0]
        base_name = restart_file.split('-1.restart')[0]

        subdirectory_name = (
            f"ph_{base_name}_{supercell_str}_{cell_tipe}_{os.path.basename(main_dir)}"
        )
        subdirectory_path = os.path.join(main_dir, subdirectory_name)
        print(f"Processing subdirectory: {subdirectory_path}")

        try:
            if os.path.exists(subdirectory_path):
                shutil.rmtree(subdirectory_path)
                print(f"Deleted existing subdirectory: {subdirectory_path}")

            os.makedirs(subdirectory_path)
            print(f"Created new subdirectory: {subdirectory_path}")

            destination_input_file = os.path.join(subdirectory_path, f"{base_name}.inp")
            shutil.copy(os.path.join(main_dir, restart_file), destination_input_file)
            print(f"Copied {restart_file} to {destination_input_file}")

            phonon_script = os.path.expanduser("~/scripts/pcp2k.py")
            log_file = "pre.log"
            command = (
                 f"python {phonon_script} -m pre -ct {cell_tipe} "
                 f"-s {' '.join(map(str, supercell_size))} -t {tolerance} > {log_file} 2>&1 &"
            )
            run_cmd(f"cd {subdirectory_path} && {command}")


        except Exception as e:
            print(f"Error processing {main_dir}: {e}")

def submit_phonopy():
    """
    Submits phonopy calculations by finding subdirectories that start with 'ph_'
    inside each main directory.
    """
    print("Submitting phonopy calculations...")
    main_directories = get_directory_list()
    if not main_directories:
        print("No main directories found for submission.")
        return

    submit_command = os.path.expanduser("~/scripts/nodes_and_send.sh")

    for main_dir in main_directories:
        try:
            subdirs = [
                name for name in os.listdir(main_dir)
                if os.path.isdir(os.path.join(main_dir, name)) and name.startswith("ph_")
            ]

            if not subdirs:
                print(f"No 'ph_' subdirectories in {main_dir}. Skipping.")
                continue

            for sub in subdirs:
                subdirectory_path = os.path.join(main_dir, sub)
                print(f"Submitting job for {subdirectory_path} using {submit_command}")
                run_cmd(f"cd {subdirectory_path} && {submit_command}")

        except Exception as e:
            print(f"Error processing {main_dir}: {e}")

def phonopy_post():
    """Executes the phonopy post-processing and plotting script in each folder."""
    print("Starting phonopy post-processing...")
    directories = get_directory_list()
    if not directories:
        print("No directories found for processing.")
        return

    initial_dir = os.getcwd()

    for dir in directories:
        subfolders = [
            os.path.join(dir, f)
            for f in os.listdir(dir)
            if f.startswith('ph_') and os.path.isdir(os.path.join(dir, f))
        ]

        if not subfolders:
            print(f"No 'ph_' subfolders found in {dir}. Skipping.")
            continue

        for subfolder in subfolders:
            # Check if FORCE_SETS file exists in the subfolder, skip if it is present
            #force_sets_path = os.path.join(subfolder, "FORCE_SETS")
            #if os.path.isfile(force_sets_path):
            #    print(f"FORCE_SETS found in {subfolder}, skipping.")
            #    continue     

            print(f"Processing: {subfolder}")
            try:
                os.chdir(subfolder)

                script_path = os.path.expanduser("~/scripts/pcp2k.py")
                log_file = "post.log"

                command_post = (
                    f"python {script_path} -m post -tp > {log_file} 2>&1"
                )
                print(f"Running: {command_post}")
                run_cmd(command_post)

                print(f"Finished post-processing and plotting in {subfolder}")

            except Exception as e:
                print(f"Error processing {subfolder}: {e}")
            finally:
                os.chdir(initial_dir)


def qha_post(temp, pressure):
    """
    Performs QHA analysis by:
      1. Executing create_ev() to generate the e-v.dat file.
      2. Copying e-v.dat from the current directory to a new folder 'thermal_properties'.
      3. Finding all thermal_properties.yaml files from directories provided by get_directory_list,
         copying them to the thermal_properties folder with renamed filenames based on their parent folder.
      4. Executing the command 'phonopy-qha -p -s --tmax={temp} e-v.dat <yaml files>' inside that folder.
    """
    print("Starting QHA analysis...")

    # Generate e-v.dat file.
    create_ev()

    directories = get_directory_list()
    thermal_properties_dir = "thermal_properties"
    os.makedirs(thermal_properties_dir, exist_ok=True)

    # Copy e-v.dat from the current directory (parent directory) to the thermal_properties directory.
    parent_dir = os.getcwd()
    ev_dat_src = os.path.join(parent_dir, "e-v.dat")
    if os.path.exists(ev_dat_src):
        shutil.copy(ev_dat_src, thermal_properties_dir)
    else:
        print("Warning: e-v.dat not found in the current directory.")

    # Copy thermal_properties.yaml files into the thermal_properties directory with renamed files.
    thermal_property_files = []
    for directory in directories:
        for root, _, files in os.walk(directory):
            for file in files:
                if file == "thermal_properties.yaml":
                    file_path = os.path.join(root, file)
                    parent_folder = os.path.basename(root)
                    new_file_name = f"thermal_properties_{parent_folder}.yaml"
                    dest_path = os.path.join(thermal_properties_dir, new_file_name)
                    shutil.copy(file_path, dest_path)
                    thermal_property_files.append(new_file_name)

    if not thermal_property_files:
        print("No thermal_properties.yaml files found. Exiting QHA analysis.")
        return

    # Build the command using the new file paths (which are relative to the thermal_properties directory).
    cmd = f"phonopy-qha -p -s --tmax={temp} e-v.dat " + " ".join(thermal_property_files)
    print(f"Executing in '{thermal_properties_dir}': {cmd}")

    # Change to the thermal_properties directory, execute the command, and then return to the original directory.
    original_dir = os.getcwd()
    try:
        os.chdir(thermal_properties_dir)
        run_cmd(cmd)
        os.system(f"python ~/scripts/plotting/plot_helmholtz.py")
        os.system(f"python ~/scripts/plotting/plot_thermal_prop.py")
    finally:
        os.chdir(original_dir)
    

def pore_calc():
    """
    Converts .restart → .xyz → .cif, runs network to get .res,
    and aggregates results into pore_data.csv.
    """
    import os

    main_csv = os.path.join(os.getcwd(), "pore_data.csv")
    with open(main_csv, "w") as outf:
        outf.write(
            "Directory,Largest_Included_Sphere,"
            "Largest_Free_Sphere,Pore_Limiting_Diameter\n"
        )

    for directory in get_directory_list():
        orig = os.getcwd()
        os.chdir(directory)

        # 1) restart → xyz
        restart = "R2SCAN-1.restart"
        run_cmd(f"python ~/scripts/file_conversion/cp2k2xyz.py {restart}")

        # 2) xyz → cif
        xyz = restart.rsplit(".", 1)[0] + ".xyz"
        run_cmd(f"python ~/scripts/file_conversion/convert_xyz2cif.py {xyz}")

        # 3) pick the .cif file
        cif_files = get_cif_files(".")
        if not cif_files:
            os.chdir(orig)
            continue
        cif = next((f for f in cif_files if f.startswith("Zn_MOF5")), cif_files[0])
        print(f"Selected .cif in {directory}: {cif}", flush=True)

        # 4) network → .res
        base = os.path.splitext(cif)[0]
        res  = base + ".res"
        run_cmd(f"network -ha -res {res} {cif}")
        if not os.path.isfile(res):
            raise FileNotFoundError(f"No .res for {cif} in {directory}")

        # 5) read & append
        with open(res) as r:
            struct, lis, lfs, pld = r.readline().split()
        lis, lfs, pld = map(float, (lis, lfs, pld))

        with open(main_csv, "a") as outf:
            outf.write(f"{directory},{lis:.5f},{lfs:.5f},{pld:.5f}\n")

        os.chdir(orig)

    print(f"All data saved to {main_csv}.")



def main():
    parser = argparse.ArgumentParser(
        description="Script with multiple modes for processing structures."
    )
    parser.add_argument(
        '-m',
        '--mode',
        required=True,
        choices=[
            'qha_pre', 'create_ev', 'pore_calc', 'apply_strain',
            'phonopy_pre', 'phonopy_post', 'qha_post', 'submit_phonopy',
            'submit_geo'
        ],
        help="Mode of operation."
    )
    parser.add_argument(
        '-v',
        '--volume',
        type=float,
        default=250,
        help='Volume increment for strain (default: 20 Å³).'
    )
    parser.add_argument(
        '-nv',
        '--number_volumes',
        type=int,
        default=7,
        help='Number of volumes to generate on each side of original (default: 5).'
    )
    parser.add_argument(
        '-t',
        '--tolerance',
        type=float,
        default=0.00005,
        help='Tolerance value for phonopy_pre (default: 0.00005).'
    )
    parser.add_argument(
        '--temp',
        default=1000,
        type=int,
        help='Temperature for QHA (default: 400 K).'
    )
    parser.add_argument(
        '--pressure',
        default=0,
        type=int,
        help='Pressure for QHA (default: 0 GPa).'
    )
    parser.add_argument(
        '--cell_tipe',
        default='conv',
        help='Cell type for phonopy_pre (default: conv).'
    )
    parser.add_argument(
        '--supercell_size',
        nargs=3,
        type=int,
        default=[1, 1, 1],
        help='Supercell size for phonopy_pre as three integers (default: 1 1 1).'
    )

    args = parser.parse_args()
    print("Arguments used:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")

    if args.mode == 'qha_pre':
        qha_pre()
    elif args.mode == 'create_ev':
        create_ev()
    elif args.mode == 'pore_calc':
        pore_calc()
    elif args.mode == 'apply_strain':
        apply_strain(
            volume_increment=args.volume,
            num_volumes=args.number_volumes
        )
    elif args.mode == 'phonopy_pre':
        phonopy_pre(
            tolerance=args.tolerance,
            cell_tipe=args.cell_tipe,
            supercell_size=args.supercell_size
        )
    elif args.mode == 'phonopy_post':
        phonopy_post()
    elif args.mode == 'qha_post':
        qha_post(args.temp, args.pressure)
    elif args.mode == 'submit_phonopy':
        submit_phonopy()
    elif args.mode == 'submit_geo':
        submit_geo()
    else:
        print(f"Error: Unknown mode {args.mode}.")
        sys.exit(1)


if __name__ == "__main__":
    main()
