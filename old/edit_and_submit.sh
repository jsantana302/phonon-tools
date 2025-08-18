#!/bin/bash

# Set your custom folder ranges here
folders="{00..19}"

# Expand the folder ranges and store them in a variable
eval folders_expanded=($folders)

# Set the control variable to decide if jobs should be sent
send_calculation=true  # Change this to 'false' if you don't want to submit jobs

# Modify EPS_SCF, KPOINTS (SCHEME MONKHORST-PACK), and lattice vectors A, B, C under CELL and CELL_REF in .inp files
for folder in "${folders_expanded[@]}"; do
    # Modify EPS_SCF
    find "$folder" -name "*.inp" -exec sed -i '/EPS_SCF/c\         EPS_SCF 1e-09' {} +

    # Modify KPOINTS
    find "$folder" -name "*.inp" -exec sed -i '/SCHEME MONKHORST-PACK/c\         SCHEME MONKHORST-PACK 4 4 4' {} +

    # Modify lattice vectors A, B, C under the &CELL block (26.074490978208484) and &CELL_REF block (31.28731439998339)
    awk '
    BEGIN {cell_found=0; cell_ref_found=0}
    /&CELL/ {cell_found=1; cell_ref_found=0}
    /&CELL_REF/ {cell_found=0; cell_ref_found=1}
    /&END CELL_REF/ {cell_ref_found=0}

    # Replace lines under &CELL
#    cell_found && /^ *A / {$0="         A 26.074490978208484 0.0 0.0"}
#    cell_found && /^ *B / {$0="         B 0.0 26.074490978208484 0.0"}
#    cell_found && /^ *C / {$0="         C 0.0 0.0 26.074490978208484"}

    # Replace lines under &CELL_REF
    cell_ref_found && /^ *A / {$0="            A 30.0 0.0 0.0"}
    cell_ref_found && /^ *B / {$0="            B 0.0 30.0 0.0"}
    cell_ref_found && /^ *C / {$0="            C 0.0 0.0 30.0"}

    {print $0}
    ' "$folder"/*.inp > temp && mv temp "$folder"/*.inp

done

# Copy the run_supercel.job file to each selected folder
for folder in "${folders_expanded[@]}"; do
    cp run_supercel.job "$folder/"
done

# Submit the jobs in each selected folder using sbatch, but only if send_calculation is true
if [ "$send_calculation" = true ]; then
    for folder in "${folders_expanded[@]}"; do
        (cd "$folder" && sbatch run_supercel.job)
    done
else
    echo "Jobs not submitted. Set send_calculation to true if you want to submit."
fi
