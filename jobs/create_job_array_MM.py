#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import argparse


def main() -> None:
    p = argparse.ArgumentParser(
        description="Write a Slurm CP2K array job for M-###.inp files."
    )
    p.add_argument("total", type=int, help="total number of input indices")
    p.add_argument("parallel", type=int, help="max concurrent array tasks")
    p.add_argument("--pad", type=int, default=3, help="zero pad width")
    p.add_argument("--prefix", default="M-", help="input prefix")
    p.add_argument("--ext", default=".inp", help="input extension")
    p.add_argument("--time", default="24:00:00", help="walltime")
    args = p.parse_args()

    if args.total < 1 or args.parallel < 1:
        p.error("total and parallel must be >= 1")
    if args.pad < 1:
        p.error("pad must be >= 1")

    tmpl = (
        "#!/bin/bash\n"
        "#SBATCH --partition=standard96s\n"
        f"#SBATCH --time={args.time}\n"
        "#SBATCH --nodes=3\n"
        "#SBATCH --ntasks-per-node=48\n"
        "#SBATCH --cpus-per-task=2\n"
        f"#SBATCH --array=1-{args.total}%{args.parallel}\n"
        "#SBATCH --requeue\n\n"
        "export PREFERRED_SOFTWARE_STACK=nhr-lmod\n"
        "source /sw/etc/profile/profile.sh\n"
        "module load cp2k/2024.1\n\n"
        "ulimit -s unlimited\n"
        "export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n"
        "export SLURM_CPU_BIND=none\n\n"
        f'idx=$(printf "%0{args.pad}d" "$SLURM_ARRAY_TASK_ID")\n'
        f'inp="{args.prefix}${{idx}}{args.ext}"\n'
        'if [[ ! -f "$inp" ]]; then\n'
        '  echo ">> $inp not found"; exit 0\n'
        "fi\n"
        f'out="${{inp%{args.ext}}}.out"\n'
        'echo "=== Running $inp ==="\n'
        'srun cp2k.psmp -i "$inp" -o "$out"\n'
        'find . -type f -name "*RESTART.kp.bak*" -delete\n'
    )

    with open("cp2k_array.job", "w", encoding="utf-8") as f:
        f.write(tmpl)
    print("Wrote cp2k_array.job")


if __name__ == "__main__":
    main()
