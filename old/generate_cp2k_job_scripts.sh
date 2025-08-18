#!/usr/bin/env bash
set -euo pipefail

# collect all two-digit folders in sorted order
dirs=( $(ls -d [0-9][0-9] | sort) )

# for each block of 5 folders, emit one cp2k_N.job
for ((i=0; i<${#dirs[@]}; i+=5)); do
  # group index = 0,1,2,...
  grp=$(( i / 5 ))
  # slice out up to 5 folders
  block=( "${dirs[@]:i:5}" )
  cat > cp2k_${grp}.job <<EOF
#!/bin/bash
#SBATCH --partition=standard96s
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
export PREFERRED_SOFTWARE_STACK=nhr-lmod
source /sw/etc/profile/profile.sh
module load cp2k/2024.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
export SLURM_CPU_BIND=none

indices="${block[*]}"
for idx in \$indices; do
  cd \$idx
  inp_file=\$(ls *.inp | head -n1)
  srun cp2k.psmp -i "\$inp_file" > "\${inp_file%.inp}.out"
  cd ..
done

find . -type f -name "*RESTART.kp.bak*" -exec rm -f {} +
EOF
done

