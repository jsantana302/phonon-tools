#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import argparse


def main() -> None:
    p = argparse.ArgumentParser(
        description="Write a Slurm CP2K array job script."
    )
    p.add_argument("total", type=int, help="total number of jobs")
    p.add_argument("parallel", type=int, help="max concurrent jobs")
    args = p.parse_args()

    if args.total < 1 or args.parallel < 1:
        p.error("total and parallel must be >= 1")

    tmpl = (
        "#!/bin/bash\n"
        "#SBATCH --partition=standard96s\n"
        "#SBATCH --time=12:00:00\n"
        "#SBATCH --nodes=3\n"
        "#SBATCH --ntasks-per-node=48\n"
        "#SBATCH --cpus-per-task=2\n"
        "#SBATCH --array=1-{total}%{parallel}\n"
        "#SBATCH --requeue\n\n"
        "export PREFERRED_SOFTWARE_STACK=nhr-lmod\n"
        "source /sw/etc/profile/profile.sh\n"
        "module load cp2k/2024.1\n\n"
        "ulimit -s unlimited\n"
        "export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n"
        "export SLURM_CPU_BIND=none\n\n"
        'i=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")\n'
        'cd "$i" || exit 1\n'
        'inp=$(ls *.inp | head -n1)\n'
        'srun cp2k.psmp -i "$inp" > "${{inp%.inp}}.out"\n'
    )

    with open("cp2k_array.job", "w", encoding="utf-8") as f:
        f.write(tmpl.format(total=args.total, parallel=args.parallel))

    print(f"Wrote cp2k_array.job for 1-{args.total}%{args.parallel}")


if __name__ == "__main__":
    main()
