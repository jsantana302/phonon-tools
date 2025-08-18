#!/bin/bash

# Default values
DEFAULT_NUM_NODES=10
DEFAULT_NTASKS_PER_NODE=48
DEFAULT_CPUS_PER_TASK=2

# Check if arguments are provided (NUM_NODES, NTASKS_PER_NODE, CPUS_PER_TASK)
if [[ $# -ge 1 ]]; then
    NUM_NODES=$1
else
    echo "No NUM_NODES specified, defaulting to $DEFAULT_NUM_NODES"
    NUM_NODES=$DEFAULT_NUM_NODES
fi

if [[ $# -ge 2 ]]; then
    NTASKS_PER_NODE=$2
else
    echo "No NTASKS_PER_NODE specified, defaulting to $DEFAULT_NTASKS_PER_NODE"
    NTASKS_PER_NODE=$DEFAULT_NTASKS_PER_NODE
fi

if [[ $# -ge 3 ]]; then
    CPUS_PER_TASK=$3
else
    echo "No CPUS_PER_TASK specified, defaulting to $DEFAULT_CPUS_PER_TASK"
    CPUS_PER_TASK=$DEFAULT_CPUS_PER_TASK
fi

# Match files with format cp2k_{number}.job
JOB_FILES=(cp2k_*.job)

for FILE in "${JOB_FILES[@]}"; do
    if [[ -f "$FILE" ]]; then
        if [[ "$FILE" == "cp2k_1.job" ]]; then
            sed -i "s/^#SBATCH --nodes=.*/#SBATCH --nodes=12/" "$FILE"
        else
            sed -i "s/^#SBATCH --nodes=.*/#SBATCH --nodes=$NUM_NODES/" "$FILE"
        fi
        sed -i "s/^#SBATCH --ntasks-per-node=.*/#SBATCH --ntasks-per-node=$NTASKS_PER_NODE/" "$FILE"
        sed -i "s/^#SBATCH --cpus-per-task=.*/#SBATCH --cpus-per-task=$CPUS_PER_TASK/" "$FILE"
    fi
done

for FILE in "${JOB_FILES[@]}"; do
    if [[ -f "$FILE" ]]; then
        sbatch "$FILE"
    fi
done

