#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import re
import argparse
from phonon_tools import get_directory_list

def check_keyword(path, pattern):
    """
    Check if the given regex pattern is found in the file at path.
    """
    try:
        with open(path, 'r', errors='ignore') as fh:
            for line in fh:
                if pattern.search(line):
                    return True
    except OSError as e:
        print(f"Error reading {path}: {e}", flush=True)
    return False


def review_geo(working_dir=None):
    dirs = get_directory_list(working_dir)
    # compile a regex to match the key phrase, ignoring extra chars or punctuation
    keyword_pattern = re.compile(r"GEOMETRY\s+OPTIMIZATION\s+COMPLETED")
    for d in dirs:
        out_files = [f for f in os.listdir(d) if f.endswith('.out')]
        if not out_files:
            print(f"No .out files in {d}", flush=True)
            continue
        for f in out_files:
            path = os.path.join(d, f)
            status = 'OK' if check_keyword(path, keyword_pattern) else 'MISSING'
            print(f"{path}: {status}", flush=True)


def main():
    parser = argparse.ArgumentParser(prog='review_script')
    parser.add_argument('-m', '--mode', choices=['geo'], required=True,
                        help='Mode of operation; currently only "geo" is supported.')
    parser.add_argument('-d', '--working-dir',
                        help='Directory in which to scan for comp_/exp_ subdirectories.')
    args = parser.parse_args()
    if args.mode == 'geo':
        review_geo(args.working_dir)
    else:
        parser.error(f"Unknown mode: {args.mode}")

if __name__ == '__main__':
    main()
