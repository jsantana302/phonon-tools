#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import numpy as np

# Replace '340.0' with the actual temperature you are working with
T = input("Enter the temperature value: ")
 

# Load the full array from the text file
parameters_filename = f'scph_parameters_T{T}.0'
parameters_filename_npy = f'scph_parameters_T{T}.0.npy'
 
full_parameters = np.loadtxt(parameters_filename)


# Save the last row in `.npy` format
np.save(parameters_filename_npy, full_parameters)

#emula al hiphive script:
loaded_parameters = np.load(parameters_filename_npy, allow_pickle=True)
parameters_start = loaded_parameters[-1,:] 


print(f'Parametros a la T={T}')
print(parameters_start)

print('Tamaño parametros total:')
print(loaded_parameters.shape)

