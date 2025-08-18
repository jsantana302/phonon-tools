#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import sys
import numpy as np
from aim2dat.strct import StructureCollection
from aim2dat.strct.strct_manipulation import scale_unit_cell

# Define the strain percentages for compression and expansion
STRAIN_VALUES = [0.04] 
vasp = False  # Set to True to generate POSCAR files

# Main script
def main():
    # Check if an input file is provided as a command-line argument
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <input_file.xyz>")
        return
    
    # Get the input file path from the first command-line argument
    input_file = sys.argv[1]
    structure_name = os.path.splitext(os.path.basename(input_file))[0]
    
    # Load the structure
    strct = StructureCollection()
    strct.append_from_file(structure_name, input_file)
    
    for strain in STRAIN_VALUES:
        # Format the strain as a string for filenames
        strain_str = f"{strain:.2f}"

        # Apply positive strain for expansion
        exp_structure = scale_unit_cell(
            strct[0],  # Use the first structure in the collection
            scaling_factors=1 + strain,  # Expand by given strain
            change_label=True
        )
        # Save the exp structure with .xyz extension
        exp_filename = f"{structure_name}_exp_{strain_str}.xyz"
        exp_structure.to_file(exp_filename)
        
        # Apply negative strain for compression
        comp_structure = scale_unit_cell(
            strct[0],  # Use the first structure in the collection
            scaling_factors=1 - strain,  # Compress by given strain
            change_label=True
        )
        # Save the comp structure with .xyz extension
        comp_filename = f"{structure_name}_comp_{strain_str}.xyz"
        comp_structure.to_file(comp_filename)

        # Convert to POSCAR format if vasp is True
        if vasp:
			#OG POSCAR
            os.system(f"python ~/scripts/get_POSCAR.py {structure_name}.xyz POSCAR")

            # Expansion POSCAR
            exp_poscar = f"POSCAR_{strain_str}_exp"
            os.system(f"python ~/scripts/get_POSCAR.py {exp_filename} {exp_poscar}")
            
            # Compression POSCAR
            comp_poscar = f"POSCAR_{strain_str}_comp"
            os.system(f"python ~/scripts/get_POSCAR.py {comp_filename} {comp_poscar}")

if __name__ == "__main__":
    main()
