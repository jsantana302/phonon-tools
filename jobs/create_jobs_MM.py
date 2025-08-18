#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import sys

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} TOTAL_INDICES")
        sys.exit(1)

    try:
        total = int(sys.argv[1])
    except ValueError:
        print("TOTAL_INDICES must be an integer")
        sys.exit(1)

    if total < 1:
        print("TOTAL_INDICES must be at least 1")
        sys.exit(1)

    max_per_job = 10
    # decide # of jobs by rounding half up: e.g. (81+5)//10 = 8
    jobs = (total + max_per_job//2) // max_per_job
    jobs = max(1, jobs)

    base, rem = divmod(total, jobs)

    template = (
        "#!/bin/bash\n"
        "#SBATCH --partition=standard96s\n"
        "#SBATCH --nodes=3\n"
        "#SBATCH --ntasks-per-node=48\n"
        "#SBATCH --cpus-per-task=2\n"
        "#SBATCH --time=24:00:00\n\n"
        "export PREFERRED_SOFTWARE_STACK=nhr-lmod\n"
        "source /sw/etc/profile/profile.sh\n"
        "module load cp2k/2024.1\n\n"
        "ulimit -s unlimited\n"
        "export OMP_NUM_THREADS=1\n"
        "export SLURM_CPU_BIND=none\n\n"
        "# <-- set the MPOSCAR indices you want here:\n"
        "indices=({indices})\n\n"
        "for i in \"${{indices[@]}}\"; do\n"
        "    idx=$(printf \"%03d\" \"$i\")\n"
        "    inp=\"M-${{idx}}.inp\"\n"
        "    if [[ ! -f $inp ]]; then\n"
        "        echo \">> $inp not found, skipping\"\n"
        "        continue\n"
        "    fi\n"
        "    echo \"=== Running $inp ===\"\n"
        "    srun cp2k.psmp -i \"$inp\" -o \"${{inp%.inp}}.out\"\n"
        "done\n\n"
        "find . -type f -name \"*RESTART.kp.bak*\" -delete\n"
    )

    start = 1
    for j in range(jobs):
        count = base + (1 if j < rem else 0)
        end = start + count - 1
        inds = " ".join(str(i) for i in range(start, end + 1))
        fname = f"cp2k_{j}.job"
        with open(fname, "w") as f:
            f.write(template.format(indices=inds))
        print(f"Wrote {fname} with indices {start}–{end}")
        start = end + 1

if __name__ == "__main__":
    main()

