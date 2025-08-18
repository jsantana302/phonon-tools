#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import shutil

# List all files in the current directory
files = os.listdir('.')

# Define the folder structure and file modification
for file in files:
    if file.endswith('.xyz') and file != 'Zn_MOF5_H_comp_0.01.xyz':
        # Extract type (comp or exp) and value from the file name
        if 'comp' in file:
            prefix = 'comp'
        elif 'exp' in file:
            prefix = 'exp'
        else:
            continue

        value = file.split('_')[-1].replace('.xyz', '')

        # Ensure folder names are prefixed properly as comp_xxx or exp_xxx
        folder_name = f"{prefix}_{value}"

        # Create the folder if it doesn't already exist
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)

        # Move the file to the corresponding folder
        shutil.move(file, os.path.join(folder_name, file))

        # Copy PBE.inp from comp_0.01 to the new folder
        source_inp = os.path.join('comp_0.01', 'PBE.inp')
        dest_inp = os.path.join(folder_name, 'PBE.inp')
        if os.path.exists(source_inp):
            shutil.copy(source_inp, dest_inp)

            # Modify COORD_FILE_NAME in the copied PBE.inp
            with open(dest_inp, 'r') as inp_file:
                lines = inp_file.readlines()

            with open(dest_inp, 'w') as inp_file:
                for line in lines:
                    if 'COORD_FILE_NAME' in line:
                        inp_file.write(f"        COORD_FILE_NAME {file}\n")
                    else:
                        inp_file.write(line)

        # Run the script in the newly created folder
        os.system(f"cd {folder_name} && ~/scripts/cell_script.bash")

print("All operations completed successfully.")
