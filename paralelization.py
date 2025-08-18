#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import re

def replace_sbatch_lines_in_files():
    # Define las líneas a reemplazar
    old_lines = [
        "#SBATCH --ntasks-per-node=96",
        "#SBATCH --cpus-per-task=1"
    ]
    new_lines = [
        "#SBATCH --ntasks-per-node=48",
        "#SBATCH --cpus-per-task=2"
    ]

    # Busca archivos que empiecen por 'cp2k_'
    for filename in os.listdir('.'):  # Recorre el directorio actual
        if filename.startswith("cp2k_") and os.path.isfile(filename):
            with open(filename, 'r') as file:
                content = file.readlines()

            # Reemplaza las líneas
            new_content = []
            for line in content:
                if old_lines[0] in line:
                    new_content.append(new_lines[0] + "\n")
                elif old_lines[1] in line:
                    new_content.append(new_lines[1] + "\n")
                else:
                    new_content.append(line)

            # Escribe el archivo con las líneas modificadas
            with open(filename, 'w') as file:
                file.writelines(new_content)

    print("Sustitución completada en los archivos 'cp2k_'.")

if __name__ == "__main__":
    replace_sbatch_lines_in_files()

