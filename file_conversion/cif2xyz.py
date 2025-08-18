#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import glob
from ase.io import read, write


def main():
    for cif in glob.glob("*.cif"):
        atoms = read(cif)
        write(f"{cif[:-4]}.xyz", atoms)


if __name__ == "__main__":
    main()

