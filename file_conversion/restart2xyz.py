#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
from ase.io import read, write

# List all files in the current directory
files = os.listdir()

# Loop through files and process .inp files
for file in files:
    if file.endswith(".restart"):
        # Create output file name by replacing .inp with .xyz
        output_xyz = file.replace(".restart", ".xyz")

        print(f"Processing {file}...")

        # Read the CP2K input file
        try:
            atoms = read(file, format="cp2k-restart")
            write(output_xyz, atoms)
            print(f"Converted {file} to {output_xyz}")
        except Exception as e:
            print(f"Failed to convert {file}: {e}")

print("Conversion completed.")

