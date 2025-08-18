#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import subprocess
import nglview as nv
from ase.io import read
import spglib
import re

def get_phonopy_file(base_path):
    phonopy_yaml = os.path.join(base_path, "phonopy.yaml")
    phonopy_disp_yaml = os.path.join(base_path, "phonopy_disp.yaml")
    
    if os.path.exists(phonopy_yaml):
        return phonopy_yaml
    elif os.path.exists(phonopy_disp_yaml):
        return phonopy_disp_yaml
    else:
        raise FileNotFoundError(f"Neither 'phonopy.yaml' nor 'phonopy_disp.yaml' found in {base_path}")

def get_force_file(base_path):
    force_sets = os.path.join(base_path, "FORCE_SETS")
    force_constants = os.path.join(base_path, "FORCE_CONSTANTS")
    
    if os.path.exists(force_sets):
        return {"force_sets_file_name": "FORCE_SETS"}
    elif os.path.exists(force_constants):
        return {"force_constants_file_name": "FORCE_CONSTANTS"}
    else:
        raise FileNotFoundError(f"Neither 'FORCE_SETS' nor 'FORCE_CONSTANTS' found in {base_path}")

#For strain and qha

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



def get_directory_list(working_dir=None):
    """Return comp_ and exp_ dirs, sorted by their numeric suffix."""
    if working_dir is None:
        working_dir = os.getcwd()

    print(f"Scanning for directories in {working_dir}...", flush=True)

    # only match names like comp_1.0 or exp_2
    pattern = re.compile(r'^(comp|exp)_(\d+(?:\.\d+)?)$')
    matched = []
    for d in os.listdir(working_dir):
        path = os.path.join(working_dir, d)
        if not os.path.isdir(path):
            continue
        m = pattern.match(d)
        if m:
            matched.append((m.group(1), float(m.group(2)), d))

    # split by type and sort
    comp = [d for t, num, d in matched if t == 'comp']
    comp.sort(key=lambda x: float(x.split('_')[1]), reverse=True)

    exp = [d for t, num, d in matched if t == 'exp']
    exp.sort(key=lambda x: float(x.split('_')[1]))

    ordered = comp + exp
    print(f"Found directories: {ordered}", flush=True)
    return ordered


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


def get_xyz_files(directory):
    """Finds and returns a sorted list of .xyz files in the given directory."""
    xyz_files = sorted([f for f in os.listdir(directory) if f.endswith('.xyz')])
    if not xyz_files:
        print(f"Error: No XYZ files found in {directory}.")
        return None
    return xyz_files

def get_restart_files(directory):
    """Finds and returns a sorted list of restart files in the given directory."""
    restart_files = sorted([f for f in os.listdir(directory) if f.endswith('-1.restart')])
    if not restart_files:
        print(f"Error: No restart files found in {directory}.")
        return None
    return restart_files

def get_inp_files(directory):
    """Finds and returns a sorted list of .inp files in the given directory."""
    inp_files = sorted([f for f in os.listdir(directory) if f.endswith('.inp')])
    if not inp_files:
        print(f"Error: No .inp files found in {directory}.")
        return None
    return inp_files

def get_cif_files(directory):
    """Finds and returns a sorted list of .cif files in the given directory."""
    cif_files = sorted([f for f in os.listdir(directory) if f.endswith('.cif')])
    if not cif_files:
        print(f"Error: No .cif files found in {directory}.")
        return None
    return cif_files


#pcp2k
def get_default_name():
    for file in os.listdir('.'):
        if file.endswith('.inp') and file not in ['Punitcell.inp', 'Bunitcell.inp']:
            return file.rsplit('.inp', 1)[0]
    return 'MOF-5'  # Default to 'MOF-5' if no suitable .inp file is found


def create_array_job(total=19, parallel=5, outfile="cp2k_array.job"):
    if total < 0 or parallel < 1:
        raise ValueError("total must be >= 0 and parallel must be >= 1")

    tmpl = (
        "#!/bin/bash\n"
        "#SBATCH --partition=standard96s\n"
        "#SBATCH --time=12:00:00\n"
        "#SBATCH --nodes=3\n"
        "#SBATCH --ntasks-per-node=48\n"
        "#SBATCH --cpus-per-task=2\n"
        "#SBATCH --array=0-{total}%{parallel}\n"
        "#SBATCH --requeue\n\n"
        "export PREFERRED_SOFTWARE_STACK=nhr-lmod\n"
        "source /sw/etc/profile/profile.sh\n"
        "module load cp2k/2024.1\n\n"
        "ulimit -s unlimited\n"
        "export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n"
        "export SLURM_CPU_BIND=none\n\n"
        'i=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")\n'
        'cd "$i" || exit 1\n'
        'inp=$(ls *.inp | head -n1)\n'
        'srun cp2k.psmp -i "$inp" > "${{inp%.inp}}.out"\n'
    )
    with open(outfile, "w") as f:
        f.write(tmpl.format(total=total, parallel=parallel))
    print(f"Wrote {outfile} for 0-{total}%{parallel}")

def create_band_pdos_conf(outfile="band-pdos.conf",
                          dim=(1, 1, 1),
                          primitive_axis="AUTO",
                          band=None,
                          band_labels=None,
                          mesh=(3, 3, 3),
                          gamma_center=True,
                          pdos="AUTO",
                          band_points=51,
                          fc_symmetry=True):
    """Create a Phonopy band-pdos.conf file.

    Parameters
    ----------
    outfile : str
        Name of the configuration file to write.
    dim : tuple of int
        Supercell dimension used for force constants.
    primitive_axis : str
        Primitive axis setting (default: 'AUTO').
    band : list of list of float
        High-symmetry path in reciprocal space. Each element is a
        list [kx, ky, kz]. If None, a default MOF-5 path is written.
    band_labels : list of str
        Labels for the band path points. If None, defaults to
        ["$\\Gamma$", "X", "W", "K", "$\\Gamma$", "L", "U"].
    mesh : tuple of int
        Mesh grid dimensions.
    gamma_center : bool
        Whether to gamma-center the mesh.
    pdos : str
        PDOS setting (e.g. 'AUTO').
    band_points : int
        Number of points along each band segment.
    fc_symmetry : bool
        Whether to symmetrize force constants.

    Returns
    -------
    None
        Writes the configuration file to `outfile`.
    """
    if band is None:
        band = [
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.5],
            [0.5, 0.25, 0.75],
            [0.375, 0.375, 0.75],
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.625, 0.25, 0.625],
        ]
    if band_labels is None:
        band_labels = ["$\\Gamma$", "X", "W", "K", "$\\Gamma$", "L", "U"]

    with open(outfile, "w") as f:
        f.write(f"DIM = {dim[0]} {dim[1]} {dim[2]}\n")
        f.write(f"PRIMITIVE_AXIS = {primitive_axis}\n")
        f.write("BAND = " + "  ".join(
            [f"{x:g} {y:g} {z:g}" for x, y, z in band]) + "\n")
        f.write("BAND_LABELS= " + " ".join(band_labels) + "\n")
        f.write(f"MESH = {mesh[0]} {mesh[1]} {mesh[2]}\n")
        f.write(f"GAMMA_CENTER = {'.TRUE.' if gamma_center else '.FALSE.'}\n")
        f.write(f"PDOS = {pdos}\n")
        f.write(f"BAND_POINTS={band_points}\n")
        f.write(f"FC_SYMMETRY = {'.TRUE.' if fc_symmetry else '.FALSE.'}\n")

def detect_input_file(user_input: str | None = None) -> str:
    """Detect the CP2K input file to use.

    Parameters
    ----------
    user_input : str or None
        File provided by the user. If None, attempt auto-detection.

    Returns
    -------
    str
        Path to the detected input file.

    Raises
    ------
    FileNotFoundError
        If no input file can be detected automatically.
    """
    if user_input:
        return user_input

    if os.path.exists("Punitcell.inp"):
        print("Auto-detected input file: Punitcell.inp")
        return "Punitcell.inp"
    elif os.path.exists("Bunitcell.inp"):
        print("Auto-detected input file: Bunitcell.inp")
        return "Bunitcell.inp"

    inp_files = [f for f in os.listdir(".") if f.endswith(".inp")]
    if len(inp_files) == 1:
        print(f"Auto-detected input file: {inp_files[0]}")
        return inp_files[0]

    raise FileNotFoundError(
        "No input file specified and unable to auto-detect. "
        "Please use -i/--input_file."
    )

#notebook


def visualize_xyz(directory, dir_base=None, tolerances=None):
    # Use default values if none are provided
    if dir_base is None:
        dir_base = '/mnt/lustre-grete/usr/u12658/'
    if tolerances is None:
        tolerances = [1e-2, 0.005, 1e-3, 1e-4, 0.0004, 1e-5, 1e-6]
    
    xyz_files = get_xyz_files(directory)
    xyz_files = [os.path.join(directory, xyz) for xyz in xyz_files]

    for xyz in xyz_files:
        try:
            file_name = os.path.relpath(xyz, start=dir_base) if dir_base and xyz.startswith(dir_base) else xyz
            print(f"\nFile: {file_name}")
            atoms = read(xyz)
            results = []
            for tol in tolerances:
                dataset = spglib.get_symmetry_dataset(
                    (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers()),
                    symprec=tol
                )
                symmetry = dataset.international if dataset else "No symmetry found"
                results.append(f"{tol}: {symmetry}")
            for i in range(0, len(results), 4):
                print("\t".join(results[i:i+4]))
            view = nv.show_ase(atoms)
            try:
                from IPython.display import display
                display(view)
            except ImportError:
                print("Visualization skipped (no IPython display available).")
        except Exception as e:
            print(f"Failed to process {file_name}: {e}")
            continue


#modify scripts
def erase_coord(inp_path):
    with open(inp_path, 'r') as f:
        lines = f.readlines()
    out = []
    skip = False
    for l in lines:
        if '&COORD' in l:
            skip = True
            continue
        if skip and '&END COORD' in l:
            skip = False
            continue
        if not skip:
            out.append(l)
    with open(inp_path, 'w') as f:
        f.writelines(out)