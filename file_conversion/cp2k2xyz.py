#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import sys
from ase.io import read, write

def convert(infile):
    base, ext = os.path.splitext(infile)
    if ext not in ('.inp', '.restart'):
        print(f"Skipped {infile}: unsupported extension")
        return

    if not os.path.isfile(infile):
        print(f"Skipped {infile}: file not found")
        return

    outfile = base + '.xyz'
    try:
        atoms = read(infile, format='cp2k-restart')
        write(outfile, atoms)
        print(f'Converted {infile} → {outfile}')
    except Exception as e:
        print(f'Failed {infile}: {e}')

def main():
    if len(sys.argv) == 2:
        convert(sys.argv[1])
    elif len(sys.argv) == 1:
        # No file passed → convert all .inp files in current folder
        all_files = [f for f in os.listdir() if f.endswith('.inp') or f.endswith('.restart')]
        if not all_files:
            print("No .inp or .restart files found in this directory.")
            return
        for f in all_files:
            convert(f)
    else:
        print(f'Usage: {sys.argv[0]} <file.inp|file.restart>')
        sys.exit(1)

if __name__ == '__main__':
    main()
