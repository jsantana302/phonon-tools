#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import shutil

# Generalized parameter and values
parameter = 'TRUST_RADIUS'
parameter_values = [0.1, 0.2, 0.3]

# File names
input_file = 'PBE.inp'
job_file = 'cp2k.job'

# Extract the first word of the parameter and convert it to lowercase
param_name = parameter.split('_')[0].lower()

# Read the content of the PBE.inp file
with open(input_file, 'r') as f:
    content = f.readlines()

# Loop over each value in parameter_values
for value in parameter_values:
    # Create a directory using the param_name and current value
    directory = f'{param_name}_{value}'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Modify the parameter line in the input file
    modified_content = []
    for line in content:
        if parameter in line:
            modified_content.append(f'       {parameter} {value}\n')  # 7 spaces before parameter
        else:
            modified_content.append(line)

    # Write the modified content to a new input file
    new_input_file = f'PBE_{value}.inp'
    with open(new_input_file, 'w') as f:
        f.writelines(modified_content)

    # Move the new input file to the corresponding directory
    shutil.move(new_input_file, os.path.join(directory, new_input_file))

    # Copy the cp2k.job file to the corresponding directory
    shutil.copy(job_file, os.path.join(directory, job_file))


