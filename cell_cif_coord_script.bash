#!/usr/bin/env python3
import sys
from ase.io import read

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.cif input.inp")
        sys.exit(1)
    cif_file, inp_file = sys.argv[1], sys.argv[2]
    atoms = read(cif_file)
    cell = atoms.get_cell()
    coords = atoms.get_positions()
    syms = atoms.get_chemical_symbols()

    # build CELL section
    a, b, c = cell
    cell_lines = ["   &CELL\n",
                  f"      A {a[0]} {a[1]} {a[2]}\n",
                  f"      B {b[0]} {b[1]} {b[2]}\n",
                  f"      C {c[0]} {c[1]} {c[2]}\n",
                  "      PERIODIC XYZ\n",
                  "      MULTIPLE_UNIT_CELL 1 1 1\n",
                  "   &END CELL\n"]

    # build COORD section
    coord_lines = ["   &COORD\n"]
    for sym, pos in zip(syms, coords):
        coord_lines.append(
            f"      {sym} {pos[0]} {pos[1]} {pos[2]}\n"
        )
    coord_lines.append("   &END COORD\n")

    # read inp and strip old blocks
    with open(inp_file) as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines)
                 if l.strip().startswith("&SUBSYS"))
    end = next(i for i, l in enumerate(lines[start:], start)
               if l.strip().startswith("&END SUBSYS"))
    new_body = []
    i = start
    while i <= end:
        l = lines[i]
        if l.strip().startswith("&CELL"):
            i += 1
            while not lines[i].strip().startswith("&END CELL"):
                i += 1
            i += 1
            continue
        if l.strip().startswith("&COORD"):
            i += 1
            while not lines[i].strip().startswith("&END COORD"):
                i += 1
            i += 1
            continue
        new_body.append(l)
        i += 1

    # insert new sections
    out = lines[:start]
    for l in new_body:
        out.append(l)
        if l.strip().startswith("&SUBSYS"):
            out.extend(cell_lines)
            out.extend(coord_lines)
    out.extend(lines[end+1:])

    # write back
    with open(inp_file, "w") as f:
        f.writelines(out)
    print(f"CELL and COORD updated in {inp_file}")

if __name__ == '__main__':
    main()

