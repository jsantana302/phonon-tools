"""
This script processes a directory structure to calculate the computational time spent on
different types of calculations using SLURM output files.

Key Features:
- Identifies directories as either optimization calculations or phonon calculations.
- Distinguishes between two categories of phonon calculations:
  - Phonon Main: SLURM files found directly in the phonon calculation directory.
  - Phonon Sub: SLURM files found within numeric subdirectories (e.g., '00', '01', ..., '18').
- Skips directories that start with 'old' or are irrelevant to the calculations.
- Outputs the results in a human-readable summary format, saved to a file named 
  'calculation_time_summary.txt'.

Usage:
- Place the script in the root directory to analyze.
- Run the script: `python retrieve_time.py`.
- The results will be saved in the same directory.

Dependencies:
- The script assumes that SLURM output files contain a line with computational time in 
  the format: "Estimated Consumption: <value> core-hours".

Output:
- A summary of optimization and phonon calculation times for each relevant directory, 
  separated into the Main and Sub categories for phonon calculations.
"""
import os
import re

def extract_computational_time(slurm_file):
    """Extracts computational time from a slurm file."""
    try:
        with open(slurm_file, 'r') as file:
            for line in file:
                match = re.search(r"Estimated Consumption: ([\d.]+) core-hours", line)
                if match:
                    return float(match.group(1))
    except Exception as e:
        print(f"Error reading file {slurm_file}: {e}")
    return 0.0

def is_phonon_calculation(directory):
    """Checks if the directory is a phonon calculation."""
    subfolders = sorted([name for name in os.listdir(directory) if name.isdigit()])
    has_19_subfolders = subfolders == [f"{i:02}" for i in range(19)]
    files = os.listdir(directory)
    has_phonopy_file = "phonopy_disp.yaml" in files
    has_force_file = "FORCE_SETS" in files or "FORCE_CONSTANTS" in files
    return has_19_subfolders or has_phonopy_file or has_force_file

def process_directory(directory):
    """Processes a directory and separates optimization and phonon times."""
    results = []
    optimization_time = 0.0
    phonon_time_main = 0.0
    phonon_time_sub = 0.0


    # Sum up optimization times for `slurm-*` files in the current directory
    try:
        files = os.listdir(directory)
    except Exception as e:
        print(f"Error accessing directory {directory}: {e}")
        return results, optimization_time, phonon_time_main, phonon_time_sub

    for file in files:
        if file.startswith("slurm-"):
            slurm_path = os.path.join(directory, file)
            time = extract_computational_time(slurm_path)
            optimization_time += time

    if optimization_time > 0:
        results.append((directory, "Optimization", optimization_time))

    # Recursively process subdirectories
    for subdir in sorted(files):
        subdir_path = os.path.join(directory, subdir)
        if os.path.isdir(subdir_path):
            if is_phonon_calculation(subdir_path):
                print(f"Identified phonon calculation: {subdir_path}")

                # Main Phonon Time
                phonon_main_time = 0.0
                for file in os.listdir(subdir_path):
                    if file.startswith("slurm-"):
                        slurm_path = os.path.join(subdir_path, file)
                        time = extract_computational_time(slurm_path)
                        phonon_main_time += time

                if phonon_main_time > 0:
                    results.append((subdir_path, "Phonon Calculation Main", phonon_main_time))
                    phonon_time_main += phonon_main_time

                # Sub Phonon Time
                phonon_sub_time = 0.0
                for subsubdir in [f"{i:02}" for i in range(19)]:
                    subsubdir_path = os.path.join(subdir_path, subsubdir)
                    if os.path.isdir(subsubdir_path):
                        for file in os.listdir(subsubdir_path):
                            if file.startswith("slurm-"):
                                slurm_path = os.path.join(subsubdir_path, file)
                                time = extract_computational_time(slurm_path)
                                phonon_sub_time += time

                if phonon_sub_time > 0:
                    results.append((subdir_path, "Phonon Calculation Sub", phonon_sub_time))
                    phonon_time_sub += phonon_sub_time

            else:
                sub_results, sub_optimization_time, sub_phonon_main_time, sub_phonon_sub_time = process_directory(subdir_path)
                results.extend(sub_results)
                optimization_time += sub_optimization_time
                phonon_time_main += sub_phonon_main_time
                phonon_time_sub += sub_phonon_sub_time

    return results, optimization_time, phonon_time_main, phonon_time_sub

def main(directory="."):
    results, optimization_total, phonon_main_total, phonon_sub_total = process_directory(directory)

    output_lines = []
    header = f"{'Path':<80} {'Category':<30} {'Total Time (core-hours)':<20}"
    output_lines.append(header)
    output_lines.append("-" * len(header))

    for path, category, total_time in results:
        output_lines.append(f"{path:<80} {category:<30} {total_time:<20.1f}")

    output_lines.append("-" * len(header))
    output_lines.append(f"{'Total Optimization Time':<80} {'':<30} {optimization_total:<20.1f}")
    output_lines.append(f"{'Total Phonon Calculation Main Time':<80} {'':<30} {phonon_main_total:<20.1f}")
    output_lines.append(f"{'Total Phonon Calculation Sub Time':<80} {'':<30} {phonon_sub_total:<20.1f}")

    output_file = os.path.join(directory, "calculation_time_summary.txt")
    with open(output_file, "w") as f:
        for line in output_lines:
            print(line)
            f.write(line + "\n")

    print(f"\nResults saved to {output_file}")

# Run the script with the current directory as default
main(".")
