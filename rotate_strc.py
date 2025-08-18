#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import numpy as np
from aim2dat.strct import Structure
from aim2dat.strct.ext_analysis import determine_molecular_fragments
from aim2dat.strct.ext_manipulation import rotate_structure
from copy import deepcopy
from phonon_tools import get_xyz_files

# Get the list of .xyz files in the current directory
xyz_files = get_xyz_files(os.getcwd())

# Define the angles for rotation
angles = [15, 45, 90]

# Iterate over each .xyz file
for input_file in xyz_files:
    # Load the MOF structure
    MOF5 = Structure.from_file(input_file)

    # Determine molecular fragments
    fragments = determine_molecular_fragments(
        MOF5,
        exclude_elements=["Zn", "O"],
        end_point_elements="Zn",  
        cn_method="atomic_radius",
        radius_type="chen_manz",
    )

    pairs = []
    rotation_vectors = []

    # Find the furthest connected atoms in each fragment
    for frag in fragments:
        max_dist = 0
        furthest_pair = None
        furthest_vec = None

        site_indices = frag.site_attributes["parent_indices"]
        positions = frag.positions

        for idx1, pos1 in zip(site_indices, positions):
            for idx2, pos2 in zip(site_indices, positions):
                if idx1 >= idx2:
                    continue

                dist = np.linalg.norm(np.array(pos1) - np.array(pos2))
                if dist > max_dist:
                    max_dist = dist
                    furthest_pair = (idx1, idx2)
                    furthest_vec = np.array(pos2) - np.array(pos1)

        if furthest_pair:
            pairs.append(furthest_pair)
            rotation_vectors.append(furthest_vec / np.linalg.norm(furthest_vec))  # Normalize

    # Get the base name of the input file to save the rotated structures
    base_name, ext = os.path.splitext(input_file)

    # Iterate over the angles and apply rotations
    for angle in angles:
        # Make a copy of the original structure for each rotation
        rot_struct = deepcopy(MOF5)

        # Apply rotation to each linker fragment
        for pair, frag, rot_vec in zip(pairs, fragments, rotation_vectors):
            rot_struct = rotate_structure(
                rot_struct,
                site_indices=frag.site_attributes["parent_indices"],
                angles=angle,  # Degrees
                vector=rot_vec,  # Normalized rotation vector
                origin=MOF5.get_positions()[pair[0]],  # Use one of the connected atoms as origin
                change_label=False,
                wrap=True  # Ensures atoms remain in unit cell
            )

        # Save the rotated structure to a new file
        rotated_file = f"{base_name}_{angle}_degrees{ext}"
        rot_struct.to_file(rotated_file)
        print(f"Rotated structure saved to: {rotated_file}")
