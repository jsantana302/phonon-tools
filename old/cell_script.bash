#!/bin/bash

<<DOCSTRING
This script updates the &CELL section in a CP2K .inp file using lattice parameters extracted
from a .xyz file in the current directory.

Steps:
1. Finds the first .xyz and .inp files in the current directory.
2. Extracts lattice parameters from the second line of the .xyz file.
3. Removes the existing &CELL block in the .inp file within the &SUBSYS block.
4. Inserts a new &CELL block with the updated parameters.

Error Handling:
- Exits if no .xyz or .inp file is found, or if lattice parameters are invalid.

DOCSTRING

# Find the .xyz file in the current directory
xyz_file=$(ls *.xyz 2>/dev/null | head -n 1)
if [[ -z "$xyz_file" ]]; then
    echo "No .xyz file found in the directory."
    exit 1
fi

# Find the .inp file in the current directory
inp_file=$(ls *.inp 2>/dev/null | head -n 1)
if [[ -z "$inp_file" ]]; then
    echo "No .inp file found in the directory."
    exit 1
fi

# Extract the lattice line from the .xyz file and remove any unwanted properties
lattice_line=$(sed -n '2p' "$xyz_file" | sed 's/Lattice=//' | sed 's/Properties=.*//g' | tr -d '"')

# Extract lattice values from the line
lattice_values=$(echo "$lattice_line" | tr -s ' ')

# Split the lattice values into A, B, and C vectors
IFS=' ' read -r a1 a2 a3 b1 b2 b3 c1 c2 c3 <<< "$lattice_values"

# Create the CELL section with indentation
cell_section="   &CELL\n      A $a1 $a2 $a3\n      B $b1 $b2 $b3\n      C $c1 $c2 $c3\n      PERIODIC XYZ\n   &END CELL"

# Remove only the CELL section within the SUBSYS block
awk -v cell="$cell_section" '
    BEGIN { inside_subsys = 0; inside_cell = 0 }
    /&SUBSYS/ { inside_subsys = 1 }
    /&END SUBSYS/ { inside_subsys = 0 }
    inside_subsys && /&CELL/ { inside_cell = 1 }
    inside_cell && /&END CELL/ { inside_cell = 0; next }
    !(inside_subsys && inside_cell) { print }
' "$inp_file" > temp.inp && mv temp.inp "$inp_file"

# Insert the new CELL section after the &SUBSYS line
awk -v cell="$cell_section" '
    /&SUBSYS/ { print; print cell; next }
    { print }
' "$inp_file" > temp.inp && mv temp.inp "$inp_file"

echo "CELL section within &SUBSYS substituted in $inp_file."

