#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import re

# Define folder names
folders = [f"{i:02d}" for i in range(19)]

# Iterate over folders
for folder in folders:
    file_path = os.path.join(folder, f"R2SCAN-supercell-{folder}.inp")
    
    if os.path.exists(file_path):
        with open(file_path, "r") as file:
            lines = file.readlines()
        
        with open(file_path, "w") as file:
            for line in lines:
                new_line = re.sub(r"PROJECT_NAME\s+.*", f"PROJECT_NAME R2SCAN-supercell-{folder}", line)
                file.write(new_line)
        
        print(f"Updated: {file_path}")
    else:
        print(f"File not found: {file_path}")

