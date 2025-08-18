#!/bin/bash

# Usage: ./convert_cp2k_to_xyz.sh input_file.inp output_file.xyz

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file.inp output_file.xyz"
    exit 1
fi

input_file="$1"
output_file="$2"

# Extract lattice parameters from &CELL section
tmp=$(grep -A 3 "&CELL$" "$input_file" | tail -n +2 | awk '{printf "%15.12f %15.12f %15.12f ", $2, $3, $4}' | sed -e 's/    ^ *//g' -e 's/ *$//g')

# Extract atomic coordinates from &COORD section
coord_section=$(awk '/&COORD/,/END COORD/' "$input_file" | grep -v "&COORD" | grep -v "END COORD")

# Count the number of atoms
atom_count=$(echo "$coord_section" | wc -l)

# Write to the XYZ file
{
    echo "$atom_count"
    echo "Lattice=\"$tmp\""
    echo "$coord_section"
} > "$output_file"

echo "Conversion complete. Output saved to $output_file"

