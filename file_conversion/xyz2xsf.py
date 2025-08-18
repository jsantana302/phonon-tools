#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
from ase import io
import sys

if len(sys.argv) != 3:
    print("Usage: python xyz_to_xsf.py input.xyz output.xsf")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the XYZ file
atoms = io.read(input_file, format="xyz")

# Write the structure to XSF format
io.write(output_file, atoms, format="xsf")

print(f"Converted {input_file} to {output_file} using ASE")

