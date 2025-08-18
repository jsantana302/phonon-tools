#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker

delta_volume_units_percentage = True

e_v_data = pd.read_csv('e-v-full.csv')

try:
    pore_data = pd.read_csv('pore_data.csv', usecols=['Largest_Included_Sphere', 'Largest_Free_Sphere', 'Pore_Limiting_Diameter'])
    pore_data['Directory'] = e_v_data['Directory']
    merged_data = pd.merge(e_v_data, pore_data, on='Directory', how='left')
    plot_pore_data = True
except FileNotFoundError:
    merged_data = e_v_data
    plot_pore_data = False

cell_volume = merged_data["Cell volume [Å³]"].values
total_energy = merged_data["Total energy [eV]"].values
directory_name = merged_data["Directory"].values

lis = lfs = pld = None

if plot_pore_data:
    lis = merged_data["Largest_Included_Sphere"].values
    lfs = merged_data["Largest_Free_Sphere"].values
    pld = merged_data["Pore_Limiting_Diameter"].values

print(merged_data)

reference_index = np.where(directory_name == "comp_0")[0][0]
reference_energy = total_energy[reference_index]
reference_volume = cell_volume[reference_index]

delta_energy_eV = total_energy - reference_energy
delta_volume = cell_volume - reference_volume

if delta_volume_units_percentage:
    delta_volume = 100 * delta_volume / reference_volume
    x_label = "ΔV [%]"
else:
    x_label = "ΔV [Å³]"

plt.rc('font', family='serif', size=16)

fig, ax1 = plt.subplots(figsize=(12, 6), dpi=300)

ax1.plot(delta_volume, delta_energy_eV, '-o', label="ΔE vs ΔV", color='purple')
ax1.set_xlabel(x_label, fontsize=18)
ax1.set_ylabel("ΔE [eV]", fontsize=18, color='purple')
ax1.tick_params(axis='y', labelcolor='purple', labelsize=16, direction='in')
ax1.tick_params(axis='x', labelsize=16, direction='in')
ax1.grid(False)

for i, label in enumerate(directory_name):
    if "comp" in label and label != "comp_0":
        ax1.annotate(label, (delta_volume[i], delta_energy_eV[i]), fontsize=12, color='purple', xytext=(-4, 8), textcoords='offset points')
    elif label == "comp_0":
        ax1.annotate("Eq.", (delta_volume[i], delta_energy_eV[i]), fontsize=12, color='purple', xytext=(-8, 6), textcoords='offset points')
    elif "exp" in label:
        ax1.annotate(label, (delta_volume[i], delta_energy_eV[i]), fontsize=12, color='purple', xytext=(-40, 8), textcoords='offset points')

if plot_pore_data:
    ax2 = ax1.twinx()
    ax2.plot(delta_volume, lis, '-x', label="Largest Included Sphere (LIS)", color='teal')
    ax2.set_ylabel("LIS [Å]", fontsize=18, color='teal')
    ax2.grid(False)
    ax2.tick_params(axis='y', labelcolor='teal', labelsize=16, direction='in')

plt.title("EOS MOF-5", fontsize=20)
fig.tight_layout()
plt.savefig("EOS.png", dpi=300)
plt.show()
