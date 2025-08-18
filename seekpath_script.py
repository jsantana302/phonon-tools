#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import ase.io
import seekpath
from ase import Atoms

# read extended XYZ; ensure your input.xyz has a comment line like:
# Lattice="a1 a2 a3 b1 b2 b3 c1 c2 c3" Properties=species:S:1:pos:R:3
structure = ase.io.read("HKUST-1_conv.xyz", format="extxyz")
cell = structure.get_cell().tolist()
pos = structure.get_scaled_positions().tolist()
nums = structure.get_atomic_numbers().tolist()
res = seekpath.get_path((cell, pos, nums))
latt = res["primitive_lattice"]
pts = res["primitive_positions"]
tys = res["primitive_types"]
prim = Atoms(cell=latt, scaled_positions=pts, numbers=tys, pbc=True)
ase.io.write("primitive.xyz", prim, format="extxyz")

