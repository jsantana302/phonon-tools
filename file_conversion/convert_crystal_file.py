#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import sys
import os
from ase.io import read, write

# Check if the filename is provided as an argument
if len(sys.argv) != 2:
    print("Usage: python convert_xsf_to_xyz.py <filename.xsf>")
    sys.exit(1)

# Get the input filename and verify it has the correct extension
input_file = sys.argv[1]
if not input_file.endswith(".xsf"):
    print("Error: The file should have a .xsf extension.")
    sys.exit(1)

# Create output filename by replacing .xsf with .xyz
output_file = os.path.splitext(input_file)[0] + ".xyz"

# Load the XSF file and write to XYZ format
atoms = read(input_file)
write(output_file, atoms)

print(f"Conversion complete: '{input_file}' -> '{output_file}'")



