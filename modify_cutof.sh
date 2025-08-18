#!/bin/bash

# Check if the script received exactly two arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <new_cutoff_value> <new_rel_cutoff_value>"
    exit 1
fi

# Assign the arguments to variables
new_cutoff=$1
new_rel_cutoff=$2

# Loop through directories 00 to 19
for dir in $(seq -w 0 19); do
    # Check if the directory exists
    if [ -d "$dir" ]; then
        # Find all .inp files in the directory
        for file in "$dir"/*.inp; do
            # Check if .inp file exists
            if [ -f "$file" ]; then
                # Use sed to replace the cutoff and rel_cutoff values
                sed -i \
                    -e "s/^\s*CUTOFF\s*[0-9.]*\s*$/         CUTOFF $new_cutoff/" \
                    -e "s/^\s*REL_CUTOFF\s*[0-9.]*\s*$/         REL_CUTOFF $new_rel_cutoff/" \
                    "$file"
                echo "Updated $file with CUTOFF=$new_cutoff and REL_CUTOFF=$new_rel_cutoff"
            fi
        done
    fi
done

