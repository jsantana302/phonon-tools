#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import subprocess
import argparse
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(
        description='Run ModeMap in background with defaults'
    )
    p.add_argument(
        '--cell', default='R2SCAN.inp',
        help='cell input file'
    )
    p.add_argument(
        '--dim', default='1 1 1',
        help='dimension string'
    )
    g = p.add_mutually_exclusive_group()
    g.add_argument(
        '--map_1d', action='store_true',
        help='generate 1D map'
    )
    g.add_argument(
        '--map_2d', action='store_true',
        help='generate 2D map'
    )
    p.add_argument(
        '--mode', default='0 0 0 1',
        help='mode for first dim'
    )
    p.add_argument(
        '--mode_2', default='0 0 0 2',
        help='mode for second dim (2D only)'
    )
    p.add_argument(
        '--q_range', default='-1.2 1.2 0.4',
        help='q_range for first dim'
    )
    p.add_argument(
        '--q_range_2', default='-1.2 1.2 0.4',
        help='q_range for second dim'
    )
    p.add_argument(
        '--supercell', default='1 1 1',
        help='supercell string'
    )
    p.add_argument(
        '--no-cp2k', dest='cp2k',
        action='store_false',
        help='disable cp2k flag'
    )
    p.set_defaults(cp2k=True)
    return p.parse_args()

def main():
    args = parse_args()
    # prepare folder & copy inputs
    Path('modemap').mkdir(exist_ok=True)
    os.chdir('modemap')
    for fname in ('FORCE_SETS', args.cell):
        src = Path('..') / fname
        subprocess.run(['cp', str(src), '.'], check=True)
    # build ModeMap command
    cmd = [
        'python',
        str(Path.home() / 'scripts' / 'ModeMap' / 'ModeMap.py'),
        '--cell', args.cell,
        '--dim', args.dim,
        '--mode', args.mode,
        '--q_range', args.q_range,
        '--supercell', args.supercell,
    ]
    # include 2D options unless user asked for 1D
    if args.map_2d or not args.map_1d:
        cmd += [
            '--map_2d',
            '--mode_2', args.mode_2,
            '--q_range_2', args.q_range_2,
        ]
    if args.cp2k:
        cmd.append('--cp2k')
    # launch in background
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        preexec_fn=os.setsid
    )
    print(f'ModeMap started in background (PID {p.pid})')

if __name__ == '__main__':
    main()
