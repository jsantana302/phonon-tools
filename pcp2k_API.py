#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import shutil
import subprocess
import argparse
from phonon_tools import get_default_name
from phonopy.interface.cp2k import read_cp2k, write_cp2k, write_cp2k_by_filename
from phonopy import Phonopy
import numpy as np
import yaml
from phonopy.structure.atoms import PhonopyAtoms
# Import the obsolete YAML writers provided by your installation:
from phonopy.file_IO import write_disp_yaml_from_dataset, write_disp_yaml



def run_cmd(cmd):
    """Run a shell command and print its output line by line.

    If the command returns a non-zero exit code, print the error.
    """
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True
    )
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
    """Create directories for supercell files, move the files, and update job
    scripts with the correct input file names.

    Args:
        supercell_files: List of supercell file names.
        cp2k_dir: Directory containing CP2K job scripts.
        mat: Material name used in naming conventions.
    """
    shutil.copy(os.path.join(cp2k_dir, 'run_supercel.job'), '.')
    os.makedirs('00', exist_ok=True)
    supercell_00_file = f"{mat}-supercell.inp"
    if supercell_00_file in supercell_files:
        shutil.move(supercell_00_file, os.path.join('00', supercell_00_file))
        shutil.copy(os.path.join(cp2k_dir, 'run_supercel.job'), '00')
        job_file_path = os.path.join('00', 'run_supercel.job')
        input_file_base = supercell_00_file.replace('.inp', '')
        run_cmd(
            f"sed -i 's/input_file=\".*\"/input_file=\"{input_file_base}\"/' "
            f"{job_file_path}"
        )
    for supercell_file in supercell_files:
        if supercell_file != supercell_00_file:
            folder_number = supercell_file.split('-')[-1].split('.')[0]
            dir_name = f"{int(folder_number):02d}"
            os.makedirs(dir_name, exist_ok=True)
            shutil.move(supercell_file, os.path.join(dir_name, supercell_file))
            shutil.copy(os.path.join(cp2k_dir, 'run_supercel.job'), dir_name)
            job_file_path = os.path.join(dir_name, 'run_supercel.job')
            input_file_base = supercell_file.replace('.inp', '')
            run_cmd(
                f"sed -i 's/input_file=\".*\"/input_file=\"{input_file_base}\"/' "
                f"{job_file_path}"
            )


def restore_admm(input_file):
    """Restore the ADMM block in the input file and ensure that the
    MIN_PAIR_LIST_RADIUS line is present in the &QS block.

    Args:
        input_file: Name of the input file to modify.
    """
    with open(input_file, 'r') as file:
        lines = file.readlines()

    admm_lines = [
        "      &AUXILIARY_DENSITY_MATRIX_METHOD\n",
        "         ADMM_TYPE ADMMS\n",
        "      &END AUXILIARY_DENSITY_MATRIX_METHOD\n"
    ]
    min_pair_line = "         MIN_PAIR_LIST_RADIUS -1\n"
    for i, line in enumerate(lines):
        if line.strip() == "&END SCF":
            insert_idx = i + 1
            lines[insert_idx:insert_idx] = admm_lines
            break

    inside_qs = False
    for i, line in enumerate(lines):
        if line.strip() == "&QS":
            inside_qs = True
        elif line.strip() == "&END QS" and inside_qs:
            if min_pair_line.strip() not in "".join(lines).strip():
                lines.insert(i, min_pair_line)
            break

    with open(input_file, 'w') as file:
        file.writelines(lines)

def pre_mode(args, input_file, cp2k_dir, mat, tol_apply, tol_detect, cell_type, no_sym,
             supercell_size, amplitude, hybrid):
    """
    Replicates the behavior of:
      phonopy --cp2k -c {input_file} -d --dim='{supercell_size[0]} {supercell_size[1]} {supercell_size[2]}' --amplitude={amplitude}
    
    This function:
      - Preprocesses the CP2K input file.
      - Reads the CP2K input via the CP2K interface.
      - Optionally performs symmetry analysis and updates the input (writing a primitive or conventional cell).
      - Creates a Phonopy object with the given supercell dimensions and generates displacements.
      - Writes YAML files (phonopy_disp.yaml and phonopy_symcells.yaml) that mimic the CLI outputs.
      - Writes each displaced supercell as a CP2K input file.
    """

    # --- Preprocess the CP2K input file ---
    run_cmd(
        f"sed -i '/&MOTION/,/&END MOTION/d; /&TOPOLOGY/,/&END TOPOLOGY/d' {input_file}"
    )
    run_cmd(
        f"sed -i 's/EXTRAPOLATION.*/EXTRAPOLATION USE_PREV_P/' {input_file}"
    )
    run_cmd("cp /user/j.santanaandreo/u12658/scripts/jobs/phonopy_jobs/cp2k_* .")
    run_cmd(f"sed -i 's/PROJECT_NAME .*/PROJECT_NAME {mat}/' {input_file}")
    run_cmd(
        f"sed -i '/&AUXILIARY_DENSITY_MATRIX_METHOD/,/&END AUXILIARY_DENSITY_MATRIX_METHOD/d' {input_file}"
    )
    run_cmd(f"sed -i '/MIN_PAIR_LIST_RADIUS -1/d' {input_file}")

    # --- Read the CP2K input file ---
    unitcell, tree = read_cp2k(input_file)
    if not hasattr(unitcell, "copy"):
        unitcell = PhonopyAtoms(*unitcell)

    # --- (Optional) Symmetry analysis ---
    if not no_sym:
        ph_sym = Phonopy(
            unitcell,
            supercell_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            is_symmetry=True,
            symprec=tol_apply 
        )
        with open(f"sym_{tol_apply}.out", "w") as f:
            f.write(f"# Symmetry detected with symprec={tol_apply}\n")
            f.write(f"Space group type: {ph_sym.symmetry.get_international_table()}\n")
            try:
                pg_ops = ph_sym.symmetry.pointgroup_operations
                pg_value = pg_ops.get("international_table", None)
                if pg_value is None:
                    raise KeyError
            except (AttributeError, KeyError):
                pg_value = ph_sym.symmetry._pointgroup
            f.write(f"Point group: {pg_value}\n")
            f.write("\n(Add details as needed)\n")
        if cell_type == "prim":
            shutil.copy(input_file, f"{mat}_inputsym")
            write_cp2k_by_filename("Punitcell.inp", ph_sym.primitive, tree)
            shutil.move("Punitcell.inp", input_file)
        elif cell_type == "conv":
            shutil.copy(input_file, f"{mat}_inputsym")
            write_cp2k_by_filename("Bunitcell.inp", ph_sym.primitive, tree)
            shutil.move("Bunitcell.inp", input_file)

    # --- Re-read the (possibly updated) CP2K input file ---
    final_unitcell, tree_final = read_cp2k(input_file)
    if not hasattr(final_unitcell, "copy"):
        final_unitcell = PhonopyAtoms(*final_unitcell)

    # --- Create a Phonopy object to generate displacements ---
    ph = Phonopy(
        final_unitcell,
        supercell_matrix=supercell_size,
        is_symmetry=True, #busca simetría en la supercelda
        symprec=tol_detect #ESTO
    )
    ph.generate_displacements(distance=amplitude)
    
    # --- Write YAML files as produced by the CLI ---
    # For the displacement file, choose the function based on the dataset format.
    if "first_atoms" in ph.dataset:
        write_disp_yaml_from_dataset(ph.dataset, ph.supercell, filename="phonopy_disp.yaml")
    else:
        write_disp_yaml(ph.dataset, ph.supercell, filename="phonopy_disp.yaml")
    
    
    # --- Write each displaced supercell as a CP2K input file ---
    for i, scell in enumerate(ph.supercells_with_displacements):
        out_name = f"{mat}-supercell-{i:02d}.inp"
        write_cp2k_by_filename(out_name, scell, tree_final)

    # --- Optionally, restore ADMM blocks and organize files ---
    if hybrid:
        supercell_files = [f for f in os.listdir('.') 
                           if f.startswith(f"{mat}-supercell") and f.endswith(".inp")]
        for supercell_file in supercell_files:
            restore_admm(supercell_file)
    supercell_files = [f for f in os.listdir('.') 
                       if f.startswith(f"{mat}-supercell") and f.endswith(".inp")]
    create_directories_and_move_files(supercell_files, cp2k_dir, mat)

def pre_just_files_mode(mat, cp2k_dir):
    """Perform operations for the 'pre_just_files' mode.

    Create directories for each supercell file and move the files.

    Args:
        mat: Material name.
        cp2k_dir: Directory containing CP2K job scripts.
    """
    supercell_files = [
        f for f in os.listdir('.') if f.startswith(f"{mat}-supercell") and
        f.endswith(".inp")
    ]
    create_directories_and_move_files(supercell_files, cp2k_dir, mat)


def post_mode(input_file, cp2k_dir, mat, thermal_properties):
    """Perform operations for the 'post' mode.

    Process force files, run phonopy analysis and plotting, and handle
    thermal properties if requested.

    Args:
        input_file: Name of the CP2K input file.
        cp2k_dir: Directory containing CP2K job scripts.
        mat: Material name.
        thermal_properties: Boolean flag for thermal properties tasks.
    """
    supercells = []
    for dir_name in sorted(os.listdir('.')):
        if dir_name == '00':
            continue
        if (os.path.isdir(dir_name) and dir_name.isdigit() and
                len(dir_name) == 2):
            for file_name in os.listdir(dir_name):
                if file_name.endswith("-forces-1_0.xyz"):
                    supercells.append(os.path.join(dir_name, file_name))
    run_cmd("phonopy --cp2k -f " + " ".join(supercells))
    if not os.path.exists('band-pdos.conf'):
        shutil.copy(os.path.join(cp2k_dir, 'band-pdos.conf'), ".")
    run_cmd("phonopy --cp2k -c " + input_file +
            " -p -s band-pdos.conf")
    run_cmd("~/scripts/plot_phonons.py")
    if thermal_properties:
        run_cmd("phonopy --cp2k -c " + input_file +
                " -t -s band-pdos.conf")


def qha_mode(temp, pressure, n_files):
    """Perform operations for the 'qha' mode.

    Run the phonopy-qha command with the provided temperature, pressure, and
    number of thermal property files.

    Args:
        temp: Maximum temperature for QHA.
        pressure: Pressure for QHA.
        n_files: Number of thermal properties files.
    """
    files_str = " ".join([
        f"thermal_properties-{i:03d}.yaml"
        for i in range(1, n_files + 1)
    ])
    cmd = (
        f"phonopy-qha --cp2k -p -s --tmax={temp} --pressure={pressure} "
        f"e-v.dat {files_str}"
    )
    print(cmd)
    run_cmd(cmd)


def submit_mode(mat, folders, copy_run):
    """Perform operations for the 'submit' mode.

    Submit job scripts to the scheduler for specified directories.

    Args:
        mat: Material name.
        folders: List of folders to submit jobs for.
        copy_run: If True, copy the job script to each directory.
    """
    if folders:
        dirs_to_submit = [folder for folder in folders if os.path.isdir(folder)]
    else:
        dirs_to_submit = [
            dir_name for dir_name in sorted(os.listdir('.'))
            if os.path.isdir(dir_name) and dir_name.isdigit() and
            len(dir_name) == 2
        ]
    if os.path.isdir('00'):
        dirs_to_submit.insert(0, '00')
    if copy_run:
        for dir_name in dirs_to_submit:
            shutil.copy('run_supercel.job', dir_name)
    for dir_name in dirs_to_submit:
        job_name = f"{dir_name}-{mat}"
        run_cmd(f"cd {dir_name} && sbatch -J {job_name} run_supercel.job")


def main():
    """Main function to parse arguments and execute the requested mode."""
    parser = argparse.ArgumentParser(
        description='Run phonopy tasks with CP2K'
    )
    parser.add_argument(
        '-m', '--mode', required=True,
        choices=['pre', 'post', 'qha', 'submit', 'pre_just_files'],
        help=('Mode of operation: pre, post, qha, submit, or '
              'pre_just_files')
    )
    parser.add_argument(
        '-s', '--supercell_size', nargs=3, default=[2, 2, 2],
        type=int, help='Supercell sizes (required for pre)'
    )
    parser.add_argument(
        '-ct', '--cell_type', default='prim', choices=['prim', 'conv'],
        help='Cell type (required for pre)'
    )
    parser.add_argument(
        '-ta', '--tolerance_apply', default=0.00001, type=float,
        help='Symmetry tolerance to apply in the symetrization of the input structure before the supercell is generated. Default: 1e-5.'
    )
    parser.add_argument(
        '-td', '--tolerance_detection', default=0.0001, type=float,
        help='Symmetry tolerance detection to determine the number of distorted supercells needed. Default: 1e-4.'
    )
    parser.add_argument(
        '-a', '--amplitude', default=0.01, type=float,
        help='Amplitude of the distorted supercell'
    )
    parser.add_argument(
        '-n', '--name', required=False, default=get_default_name(),
        help='Material name (default: name before .inp file)'
    )
    parser.add_argument(
        '--n_files', default=9, type=int,
        help='Number of thermal properties files for QHA (default: 9)'
    )
    parser.add_argument(
        '--temp', default=1000, type=int,
        help='Temperature for QHA (default: 980)'
    )
    parser.add_argument(
        '--pressure', default=0, type=int,
        help='Pressure for QHA (default: 0)'
    )
    parser.add_argument(
        '--folders', nargs='*',
        help=('Specific folders to submit jobs for (default: all folders)')
    )
    parser.add_argument(
        '-cr', '--copy_run', action='store_true', default=True,
        help=('Copy run_supercel.job to each directory before submitting '
              '(default: true)')
    )
    parser.add_argument(
        '--no_sym', action='store_true',
        help=('Skip symmetry operations and '
              'Punitcell/Bunitcell file handling in pre mode')
    )
    parser.add_argument(
        '-tp', '--thermal_properties', action='store_true', default=False,
        help='Execute specific command in post mode (default: false)'
    )
    args = parser.parse_args()

    print(f"Mode: {args.mode}")
    print("Arguments:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")

    task = args.mode
    mat = args.name
    tol_apply = args.tolerance_apply
    tol_detect = args.tolerance_detection
    cell_type = args.cell_type
    no_sym = args.no_sym
    supercell_size = args.supercell_size
    amplitude = args.amplitude

    input_file = f"{mat}.inp"
    cp2k_dir = "/user/j.santanaandreo/u12658/files_CP2K/MOF-5"

    with open(input_file, 'r') as f:
        input_content = f.read()
    hybrid = "&AUXILIARY_DENSITY_MATRIX_METHOD" in input_content

    if task == "pre":
        pre_mode(args, input_file, cp2k_dir, mat, tol_apply, tol_detect, cell_type, no_sym,
                 supercell_size, amplitude, hybrid)
    elif task == "pre_just_files":
        pre_just_files_mode(mat, cp2k_dir)
    elif task == "post":
        post_mode(input_file, cp2k_dir, mat, args.thermal_properties)
    elif task == "qha":
        qha_mode(args.temp, args.pressure, args.n_files)
    elif task == "submit":
        submit_mode(mat, args.folders, args.copy_run)
    else:
        print("I don't understand")


if __name__ == '__main__':
    main()
