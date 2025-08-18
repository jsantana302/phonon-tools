#!/bin/bash

<<DOCSTRING
This script updates the &CELL and &COORD sections in a CP2K .inp file using lattice
parameters and atomic coordinates extracted from an .xyz file in the current directory.

Steps:
1. Finds the first .xyz and .inp files in the current directory.
2. Extracts lattice parameters from the second line of the .xyz file.
3. Extracts atomic coordinates from the .xyz file (starting from the third line).
4. Removes the existing &CELL and &COORD blocks in the .inp file (within the &SUBSYS block).
5. Inserts new &CELL and &COORD blocks with the updated parameters and coordinates.

Exits if no .xyz or .inp file is found, or if lattice parameters are invalid.
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

# Extract the lattice line (2nd line) from the .xyz file,
# remove the "Lattice=" part and any trailing properties, and squeeze spaces.
lattice_line=$(sed -n '2p' "$xyz_file" | sed 's/Lattice=//' | sed 's/Properties=.*//g' | tr -d '"')
lattice_values=$(echo "$lattice_line" | tr -s ' ')

# Read the lattice values into individual variables:
IFS=' ' read -r a1 a2 a3 b1 b2 b3 c1 c2 c3 <<< "$lattice_values"
if [[ -z "$a1" || -z "$b1" || -z "$c1" ]]; then
    echo "Error: Lattice parameters not found or invalid in $xyz_file."
    exit 1
fi

# Construct the new &CELL block. (Adjust indentation as needed.)
cell_section="   &CELL
      A $a1 $a2 $a3
      B $b1 $b2 $b3
      C $c1 $c2 $c3
      PERIODIC XYZ
      MULTIPLE_UNIT_CELL 1 1 1
   &END CELL"

# Extract the coordinate lines (lines 3 onward) from the .xyz file,
# and indent each line (adding six spaces).
coords=$(tail -n +3 "$xyz_file" | sed 's/^/      /')
coord_section="   &COORD
$coords
   &END COORD"

# Remove the existing &CELL and &COORD blocks (only within the &SUBSYS block)
awk '
    BEGIN {inside_subsys=0; inside_cell=0; inside_coord=0}
    /&SUBSYS/ {inside_subsys=1}
    /&END SUBSYS/ {inside_subsys=0}
    inside_subsys && /&CELL/ {inside_cell=1; next}
    inside_cell && /&END CELL/ {inside_cell=0; next}
    inside_subsys && /&COORD/ {inside_coord=1; next}
    inside_coord && /&END COORD/ {inside_coord=0; next}
    !(inside_cell || inside_coord) {print}
' "$inp_file" > temp.inp && mv temp.inp "$inp_file"

# Insert the new &CELL and &COORD blocks after the &SUBSYS line
awk -v cell="$cell_section" -v coord="$coord_section" '
    /&SUBSYS/ {print; print cell; print coord; next}
    {print}
' "$inp_file" > temp.inp && mv temp.inp "$inp_file"

echo "CELL and COORD sections within &SUBSYS substituted in $inp_file."
