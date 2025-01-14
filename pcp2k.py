#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python

import os
import shutil
import subprocess
import argparse

def run_cmd(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    if rc != 0:
        print(f"Error executing command: {cmd}\n")
        err = process.stderr.read()
        if err:
            print(err.strip())

def create_directories_and_move_files(supercell_files, cp2k_dir, mat):
    shutil.copy(os.path.join(cp2k_dir, 'run_supercel.job'), '.')

    # Create folder named 00 and move {mat}-supercell.inp to this folder
    os.makedirs('00', exist_ok=True)
    supercell_00_file = f"{mat}-supercell.inp"
    if supercell_00_file in supercell_files:
        shutil.move(supercell_00_file, os.path.join('00', supercell_00_file))
        shutil.copy(os.path.join(cp2k_dir, 'run_supercel.job'), '00')
        job_file_path = os.path.join('00', 'run_supercel.job')
        input_file_base = supercell_00_file.replace('.inp', '')
        run_cmd(f"sed -i 's/input_file=\".*\"/input_file=\"{input_file_base}\"/' {job_file_path}")

    for supercell_file in supercell_files:
        if supercell_file != supercell_00_file:
            # Extract the number from the supercell file name
            folder_number = supercell_file.split('-')[-1].split('.')[0]
            dir_name = f"{int(folder_number):02d}"  # Convert the number to two digits
            os.makedirs(dir_name, exist_ok=True)
            shutil.move(supercell_file, os.path.join(dir_name, supercell_file))
            shutil.copy(os.path.join(cp2k_dir, 'run_supercel.job'), dir_name)

            # Update the input file name in the job script directly
            job_file_path = os.path.join(dir_name, 'run_supercel.job')
            input_file_base = supercell_file.replace('.inp', '')
            run_cmd(f"sed -i 's/input_file=\".*\"/input_file=\"{input_file_base}\"/' {job_file_path}")

def get_default_name():
    for file in os.listdir('.'):
        if file.endswith('.inp') and file not in ['Punitcell.inp', 'Bunitcell.inp']:
            return file.rsplit('.inp', 1)[0]
    return 'MOF-5'  # Default to 'MOF-5' if no suitable .inp file is found

def restore_admm(input_file):
    """
    Restores the ADMM block after the `&END SCF` line in the input file
    and ensures `MIN_PAIR_LIST_RADIUS -1` is included in the `&QS` block.
    """
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Define the ADMM block
    admm_lines = [
        "      &AUXILIARY_DENSITY_MATRIX_METHOD\n",
        "         ADMM_TYPE ADMMS\n",
        "      &END AUXILIARY_DENSITY_MATRIX_METHOD\n"
    ]

    # Define the `MIN_PAIR_LIST_RADIUS` line
    min_pair_list_radius_line = "         MIN_PAIR_LIST_RADIUS -1\n"

    # Find the `&END SCF` block to add the ADMM block after it
    for i, line in enumerate(lines):
        if line.strip() == "&END SCF":
            insert_idx = i + 1
            lines[insert_idx:insert_idx] = admm_lines
            break

    # Ensure the `MIN_PAIR_LIST_RADIUS` line exists in the `&QS` block
    inside_qs_block = False
    for i, line in enumerate(lines):
        if line.strip() == "&QS":
            inside_qs_block = True
        elif line.strip() == "&END QS" and inside_qs_block:
            # Before ending `&QS`, ensure `MIN_PAIR_LIST_RADIUS` is present
            if min_pair_list_radius_line.strip() not in "".join(lines).strip():
                lines.insert(i, min_pair_list_radius_line)
            break

    # Write the modified lines back to the file
    with open(input_file, 'w') as file:
        file.writelines(lines)


# Argument parsing
parser = argparse.ArgumentParser(description='Run phonopy tasks with CP2K')
parser.add_argument('-m', '--mode', required=True, choices=['pre', 'post', 'qha', 'submit', 'pre_just_files'], 
                    help='Mode of operation: pre, post, qha, submit, or pre_just_files')
parser.add_argument('-s', '--supercell_size', nargs=3, default=[2, 2, 2], type=int, 
                    help='Supercell sizes (required for pre)')
parser.add_argument('-ct', '--cell_type', default='prim', choices=['prim', 'conv'], 
                    help='Cell type (required for pre)')
parser.add_argument('-t', '--tolerance', default=0.00001, type=float, help='Symmetry tolerance')
parser.add_argument('-a', '--amplitude', default=0.01, type=float, help='Amplitude of the distorted supercell')
parser.add_argument('-n', '--name', required=False, default=get_default_name(), 
                    help='Material name (default: name before .inp file)')
parser.add_argument('--n_files', default=9, type=int, help='Number of thermal properties files for QHA (default: 9)')
parser.add_argument('--temp', default=1000, type=int, help='Temperature for QHA (default: 980)')
parser.add_argument('--pressure', default=0, type=int, help='Pressure for QHA (default: 0)')
parser.add_argument('--folders', nargs='*', help='Specific folders to submit jobs for (default: all folders)')
parser.add_argument('-cr', '--copy_run', action='store_true', default=True, 
                    help='Copy run_supercel.job to each directory before submitting (default: true)')
parser.add_argument('--no_sym', action='store_true', help='Skip symmetry operations and Punitcell/Bunitcell file handling in pre mode')
parser.add_argument('-tp', '--thermal_properties', action='store_true', default=False, 
                    help='Execute specific command in post mode (default: false)')

args = parser.parse_args()

# Print the mode and arguments
print(f"Mode: {args.mode}")
print("Arguments:")
for arg, value in vars(args).items():
    print(f"  {arg}: {value}")

# Variables
task = args.mode
mat = args.name
tolerance = args.tolerance
cell_type = args.cell_type
no_sym = args.no_sym
thermal_properties = args.thermal_properties
supercell_size = args.supercell_size
amplitude = args.amplitude

input_file = f"{mat}.inp"
cp2k_dir = f"/user/j.santanaandreo/u12658/files_CP2K/MOF-5"

# Check for AUXILIARY_DENSITY_MATRIX_METHOD and set hybrid
with open(input_file, 'r') as f:
    input_content = f.read()
hybrid = "&AUXILIARY_DENSITY_MATRIX_METHOD" in input_content

if task == "pre":
    # Remove &MOTION and &TOPOLOGY sections and modify extrapolation
    run_cmd(f"sed -i '/&MOTION/,/&END MOTION/d; /&TOPOLOGY/,/&END TOPOLOGY/d' {input_file}")
    run_cmd(f"sed -i 's/EXTRAPOLATION.*/EXTRAPOLATION USE_PREV_P/' {input_file}")
    run_cmd(f"cp /user/j.santanaandreo/u12658/files_CP2K/phonopy_jobs/cp2k_* .")
    run_cmd(f"sed -i 's/PROJECT_NAME .*/PROJECT_NAME {mat}/' {input_file}")
    run_cmd(f"sed -i '/&AUXILIARY_DENSITY_MATRIX_METHOD/,/&END AUXILIARY_DENSITY_MATRIX_METHOD/d' {input_file}")
    run_cmd(f"sed -i '/MIN_PAIR_LIST_RADIUS -1/d' {input_file}")

    # Handle symmetry operations if not skipped
    if not no_sym:
        run_cmd(f"phonopy --cp2k -c {input_file} --symmetry --tolerance {tolerance} > sym_{tolerance}.out")

        if cell_type == "prim":
            shutil.copy(input_file, f"{mat}_inputsym")
            shutil.move("Punitcell.inp", input_file)
            # os.remove("Bunitcell.inp")  # Uncomment if needed
        elif cell_type == "conv":
            shutil.copy(input_file, f"{mat}_inputsym")
            shutil.move("Bunitcell.inp", input_file)
            # os.remove("Punitcell.inp")  # Uncomment if needed

    # Generate distorted supercells
    run_cmd(f"phonopy --cp2k -c {input_file} -d --dim='{supercell_size[0]} {supercell_size[1]} {supercell_size[2]}' --amplitude={amplitude}")

    # Restore ADMM lines in all supercell files if hybrid flag is enabled
    if hybrid:
        run_cmd(f"cp /user/j.santanaandreo/u12658/files_CP2K/phonopy_jobs/hybrid_jobs/cp2k_* .")
        supercell_files = [f for f in os.listdir('.') if f.startswith(f"{mat}-supercell") and f.endswith(".inp")]
        for supercell_file in supercell_files:
            restore_admm(supercell_file)
			
    # Create directories for supercell files
    supercell_files = [f for f in os.listdir('.') if f.startswith(f"{mat}-supercell") and f.endswith(".inp")]
    create_directories_and_move_files(supercell_files, cp2k_dir, mat)



# Perform operations for pre_just_files mode
elif task == "pre_just_files":
    # Create directories for each supercell and move files
    supercell_files = [f for f in os.listdir() if f.startswith(f"{mat}-supercell") and f.endswith(".inp")]
    create_directories_and_move_files(supercell_files, cp2k_dir, mat)

elif task == "post":
    supercells = []  # paths to all force files
    for dir_name in sorted(os.listdir('.')):
        if dir_name == '00':
            continue
        if os.path.isdir(dir_name) and dir_name.isdigit() and len(dir_name) == 2:
            for file_name in os.listdir(dir_name):
                if file_name.endswith("-forces-1_0.xyz"):
                    supercells.append(os.path.join(dir_name, file_name))
    # Run phonopy with the collected force files
    run_cmd("phonopy --cp2k -f " + " ".join(supercells))

    # Check if band-pdos.conf exists, copy it if not
    if not os.path.exists('band-pdos.conf'):
        shutil.copy(os.path.join(cp2k_dir, 'band-pdos.conf'), ".")

    # Run phonopy for band structure and PDOS, must be runned to generate the phonopy.yaml
    run_cmd("phonopy --cp2k -c " + input_file + " -p -s band-pdos.conf")

    run_cmd("~/scripts/plot_phonons.py")

    # Check the thermal_properties and run the additional command if it is set
    if thermal_properties:
        run_cmd("phonopy --cp2k -c " + input_file + " -t -s band-pdos.conf")
	
# Perform operations for QHA
elif task == "qha":
    cmd = f"phonopy-qha --cp2k -p -s --tmax={temp} --pressure={pressure} e-v.dat " + \
          " ".join([f"thermal_properties-{i:03d}.yaml" for i in range(1, n_files + 1)])
    print(cmd)
    run_cmd(cmd)

# Submit jobs
elif task == "submit":
    if folders:
        dirs_to_submit = [folder for folder in folders if os.path.isdir(folder)]
    else:
        dirs_to_submit = [dir_name for dir_name in sorted(os.listdir('.')) if os.path.isdir(dir_name) and dir_name.isdigit() and len(dir_name) == 2]

    # Include the 00 directory if it exists
    if os.path.isdir('00'):
        dirs_to_submit.insert(0, '00')

    if copy_run:
        for dir_name in dirs_to_submit:
            shutil.copy('run_supercel.job', dir_name)
    
    for dir_name in dirs_to_submit:
        job_name = f"{dir_name}-{mat}"
        run_cmd(f"cd {dir_name} && sbatch -J {job_name} run_supercel.job")

else:
    print("I don't understand")
