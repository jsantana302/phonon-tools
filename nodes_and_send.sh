#!/bin/bash

DEFAULT_NUM_NODES=3
DEFAULT_NTASKS_PER_NODE=48
DEFAULT_CPUS_PER_TASK=2
DEFAULT_TIME_HOURS=12

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <NUM_NODES> [<NTASKS_PER_NODE> <CPUS_PER_TASK> <TIME_HOURS>]"
    echo "No NUM_NODES, defaulting to $DEFAULT_NUM_NODES"
    NUM_NODES=$DEFAULT_NUM_NODES
else
    NUM_NODES=$1
fi

if [[ $# -lt 2 ]]; then
    echo "No NTASKS_PER_NODE, defaulting to $DEFAULT_NTASKS_PER_NODE"
    NTASKS_PER_NODE=$DEFAULT_NTASKS_PER_NODE
else
    NTASKS_PER_NODE=$2
fi

if [[ $# -lt 3 ]]; then
    echo "No CPUS_PER_TASK, defaulting to $DEFAULT_CPUS_PER_TASK"
    CPUS_PER_TASK=$DEFAULT_CPUS_PER_TASK
else
    CPUS_PER_TASK=$3
fi

if [[ $# -lt 4 ]]; then
    echo "No TIME_HOURS, defaulting to $DEFAULT_TIME_HOURS"
    TIME_HOURS=$DEFAULT_TIME_HOURS
else
    TIME_HOURS=$4
fi

JOB_FILES=(cp2k_*.job)

for FILE in "${JOB_FILES[@]}"; do
    [[ -f $FILE ]] || continue
    sed -i "s/^#SBATCH --nodes=.*/#SBATCH --nodes=$NUM_NODES/" \
        "$FILE"
    sed -i "s/^#SBATCH --ntasks-per-node=.*/#SBATCH --ntasks-per-node=$NTASKS_PER_NODE/" \
        "$FILE"
    sed -i "s/^#SBATCH --cpus-per-task=.*/#SBATCH --cpus-per-task=$CPUS_PER_TASK/" \
        "$FILE"
    sed -i "s/^#SBATCH --time=.*/#SBATCH --time=${TIME_HOURS}:00:00/" \
        "$FILE"
done

for FILE in "${JOB_FILES[@]}"; do
    [[ -f $FILE ]] && sbatch "$FILE"
done
