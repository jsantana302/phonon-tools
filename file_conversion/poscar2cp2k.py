#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import argparse
from pathlib import Path
from ase.io import read

def poscar2cp2k(input_file, template, output_file):
    structure = read(input_file, format='vasp')
    tpl = Path(template).read_text().splitlines()
    start = next(i for i, l in enumerate(tpl)
                 if l.strip().upper().startswith('&COORD'))
    end = next(i for i, l in enumerate(tpl[start+1:], start+1)
               if l.strip().upper().startswith('&END COORD'))
    pre = tpl[:start]
    post = tpl[end+1:]
    coords = ['&COORD']
    for atm in structure:
        x, y, z = atm.position
        coords.append(f'  {atm.symbol} {x:.6f} {y:.6f} {z:.6f}')
    coords.append('&END COORD')
    out_lines = pre + coords + post
    Path(output_file).write_text('\n'.join(out_lines) + '\n')

def cp2k2poscar(input_file, output_file):
    lines = Path(input_file).read_text().splitlines()
    in_cell = in_coord = False
    cell = []
    syms = []
    poss = []
    for ln in lines:
        s = ln.strip()
        up = s.upper()
        if up.startswith('&CELL'):
            in_cell = True
            continue
        if in_cell:
            if up.startswith('&END CELL'):
                in_cell = False
                continue
            if s.startswith(('A ', 'B ', 'C ')):
                parts = s.split()
                cell.append([float(p) for p in parts[1:4]])
            continue
        if up.startswith('&COORD'):
            in_coord = True
            continue
        if up.startswith('&END COORD'):
            in_coord = False
            continue
        if in_coord:
            parts = s.split()
            if len(parts) < 4:
                continue
            syms.append(parts[0])
            poss.append([float(x) for x in parts[1:4]])
    if len(cell) != 3:
        raise RuntimeError('Need 3 vectors in &CELL block')
    uniq = []
    counts = []
    for sym in syms:
        if sym in uniq:
            counts[uniq.index(sym)] += 1
        else:
            uniq.append(sym)
            counts.append(1)
    out = ['Converted from CP2K', '1.0']
    for v in cell:
        out.append('  ' + ' '.join(f'{x:.16f}' for x in v))
    out.append(' '.join(uniq))
    out.append(' '.join(str(c) for c in counts))
    out.append('Cartesian')
    for pos in poss:
        out.append('  ' + ' '.join(f'{x:.16f}' for x in pos))
    Path(output_file).write_text('\n'.join(out) + '\n')

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=['poscar2cp2k', 'cp2k2poscar'],
                   required=False, default='poscar2cp2k')
    p.add_argument('-i', '--input', dest='inp', required=False,
                   default='.',
                   help='input file or directory')
    p.add_argument('-o', '--output', dest='out', required=False, 
                   default='M-{num}.inp',
                   help='output name or pattern using "{num}"')
    p.add_argument('-t', '--template', default='R2SCAN.inp',
                   help='CP2K template for poscar2cp2k')
    p.add_argument('--pattern', default='MPOSCAR-*',
                   help='glob for batch poscar2cp2k')
    args = p.parse_args()
    if args.mode == 'poscar2cp2k':
        if not args.template:
            p.error('--template is required for poscar2cp2k')
        inp = Path(args.inp)
        if inp.is_dir():
            for poscar in sorted(inp.glob(args.pattern)):
                stem = poscar.stem
                num = stem.split('-', 1)[1]
                name = args.out.format(num=num)
                out_file = poscar.parent / name
                poscar2cp2k(str(poscar), args.template, str(out_file))
        else:
            poscar2cp2k(args.inp, args.template, args.out)
    else:
        cp2k2poscar(args.inp, args.out)

if __name__ == '__main__':
    main()
