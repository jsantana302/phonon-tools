#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import sys
from ase.io import read, write

def convert(fpath):
    out = fpath.rsplit('.xyz', 1)[0] + '.cif'
    atoms = read(fpath)
    write(out, atoms)
    print(f'Converted {fpath} to {out}')

def main():
    files = sys.argv[1:] or os.listdir()
    for f in files:
        if f.endswith('.xyz') and os.path.isfile(f):
            convert(f)
    print('Conversion completed')

if __name__ == '__main__':
    main()
