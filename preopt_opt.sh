#!/bin/bash

# Find the .inp files
file=$(ls *.inp 2>/dev/null)

# Check if no .inp file is found
if [ -z "$file" ]; then
    echo "No .inp file found in the directory."
    exit 1
fi

# Check the mode passed as argument and apply changes accordingly
if [ "$1" == "preop" ]; then
    # Apply modifications for 'preop'
#    sed -i -e "s/TZV2P/DZVP/g" \
#           -e "s/tzv2p/dzvp/g" \
    sed -i -e "s/DZVP/TZV2P/g" \
           -e "s/dzvp/tzv2p/g" \
           -e "/EPS_SCF/c\        EPS_SCF 1e-06" \
           -e "s/MAX_FORCE.*/MAX_FORCE 1.9E-05/" \
           -e "s/R2SCAN-RESTART.wfn/R2SCAN-RESTART.kp/" \
           -e "s/PBE-RESTART.wfn/PBE-RESTART.kp/" "$file"

    # Check and add KEEP_ANGLES if not present
    if ! grep -q "KEEP_ANGLES" "$file"; then
        sed -i '/TYPE DIRECT_CELL_OPT/a\    KEEP_ANGLES' "$file"
    fi

elif [ "$1" == "opt" ]; then
    # Apply modifications for 'opt'
    sed -i -e "s/DZVP/TZV2P/g" \
           -e "s/dzvp/tzv2p/g" \
           -e "/EPS_SCF/c\         EPS_SCF 1e-08" \
           -e "s/MAX_FORCE.*/MAX_FORCE 1.9E-07/" "$file"
else
    echo "Invalid mode: $1. Use 'preop' or 'opt'."
    exit 1
fi

# Inform the user of successful completion
echo "Operation $1 completed on $file."

