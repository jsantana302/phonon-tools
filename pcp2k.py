#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import shutil
import subprocess
import argparse
from pathlib import Path
import re
from phonon_tools import get_default_name


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


def tidy_supercells(pattern="*supercell*.inp") -> None:
    p = Path(".")
    for f in sorted(p.glob(pattern)):
        m = re.search(r"supercell-(\d+)\.inp$", f.name)
        d = p / (f"{int(m.group(1)):02d}" if m else "00")
        d.mkdir(exist_ok=True)
        f.rename(d / f.name)
        
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

def pre_mode(input_file, mat, tolerance, cell_type, no_sym,
             supercell_size, amplitude, hybrid):
    """Perform operations for the 'pre' mode using the CP2K interface and Phonopy API.

    This includes cleaning the input file, handling symmetry, generating
    distorted supercells, and restoring ADMM blocks if needed.

    Args:
        input_file: Name of the CP2K input file.
        mat: Material name.
        tolerance: Symmetry tolerance.
        cell_type: 'prim' or 'conv' for the cell type.
        no_sym: If True, skip symmetry operations.
        supercell_size: List of three integers defining the supercell.
        amplitude: Amplitude of the supercell distortion.
        hybrid: Boolean flag indicating if hybrid ADMM settings are used.
    """
    run_cmd(
        f"sed -i '/&MOTION/,/&END MOTION/d; /&TOPOLOGY/,/"
        f"&END TOPOLOGY/d' {input_file}"
    )
    run_cmd(
        f"sed -i 's/EXTRAPOLATION.*/EXTRAPOLATION USE_PREV_P/' "
        f"{input_file}"
    )
    run_cmd(
        "cp /user/j.santanaandreo/u12658/scripts/jobs/phonopy_jobs/"
        "cp2k_* ."
    )
    run_cmd(f"sed -i 's/PROJECT_NAME .*/PROJECT_NAME {mat}/' {input_file}")
    run_cmd(
        f"sed -i '/&AUXILIARY_DENSITY_MATRIX_METHOD/,/"
        f"&END AUXILIARY_DENSITY_MATRIX_METHOD/d' {input_file}"
    )
    run_cmd(f"sed -i '/MIN_PAIR_LIST_RADIUS -1/d' {input_file}")

    if not no_sym:
        run_cmd(
            f"phonopy --cp2k -c {input_file} --symmetry --tolerance "
            f"{tolerance} > sym_{tolerance}.out"
        )
        if cell_type == "prim":
            shutil.copy(input_file, f"{mat}_inputsym")
            shutil.move("Punitcell.inp", input_file)
        elif cell_type == "conv":
            shutil.copy(input_file, f"{mat}_inputsym")
            shutil.move("Bunitcell.inp", input_file)

    run_cmd(
        f"phonopy --cp2k -c {input_file} -d --dim='{supercell_size[0]} "
        f"{supercell_size[1]} {supercell_size[2]}' --amplitude={amplitude}"
    )

    if hybrid:
        run_cmd(
            "cp /user/j.santanaandreo/u12658/files_CP2K/phonopy_jobs/"
            "hybrid_jobs/cp2k_* ."
        )
        supercell_files = [
            f for f in os.listdir('.') if f.startswith(f"{mat}-supercell") and
            f.endswith(".inp")
        ]
        for supercell_file in supercell_files:
            restore_admm(supercell_file)

    supercell_files = [
        f for f in os.listdir('.') if f.startswith(f"{mat}-supercell") and
        f.endswith(".inp")
    ]
    tidy_supercells(pattern=f"{mat}-supercell*.inp")
    run_cmd("~/scripts/file_conversion/cp2k2xyz.py")


def pre_just_files_mode(mat):
    """Create directories for each supercell file and move the files."""
    tidy_supercells(pattern=f"{mat}-supercell*.inp")


def post_mode(input_file, cp2k_dir, thermal_properties):
    """Perform operations for the 'post' mode.

    Process force files, run phonopy analysis and plotting, and handle
    thermal properties if requested.

    Args:
        input_file: Name of the CP2K input file.
        cp2k_dir: Directory containing CP2K job scripts.
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
    run_cmd("python ~/scripts/plot_phonons_simple_dos.py")
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

def main():
    """Main function to parse arguments and execute the requested mode."""
    parser = argparse.ArgumentParser(
        description='Run phonopy tasks with CP2K'
    )
    parser.add_argument(
        '-m', '--mode', required=True,
        choices=['pre', 'post', 'qha' , 'pre_just_files'],
        help=('Mode of operation: pre, post, qha, or '
              'pre_just_files')
    )
    parser.add_argument(
        '-s', '--supercell_size', nargs=3, default=[1, 1, 1],
        type=int, help='Supercell sizes (required for pre)'
    )
    parser.add_argument(
        '-ct', '--cell_type', default='conv', choices=['prim', 'conv'],
        help='Cell type (required for pre)'
    )
    parser.add_argument(
        '-t', '--tolerance', default=0.00001, type=float,
        help='Symmetry tolerance'
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
    tolerance = args.tolerance
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
        pre_mode(input_file, mat, tolerance, cell_type, no_sym,
                 supercell_size, amplitude, hybrid)
    elif task == "pre_just_files":
        pre_just_files_mode(mat)
    elif task == "post":
        post_mode(input_file, cp2k_dir, args.thermal_properties)
    elif task == "qha":
        qha_mode(args.temp, args.pressure, args.n_files)
    else:
        print("I don't understand")


if __name__ == '__main__':
    main()
