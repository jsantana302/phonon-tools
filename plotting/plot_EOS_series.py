#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from aim2dat.plots import SimplePlot

delta_volume_units_percentage = True

# read base data
e_v_data = pd.read_csv('e-v-full.csv')
dirs = e_v_data['Directory'].values
vols = e_v_data['Cell volume [Å³]'].values
ener = e_v_data['Total energy [eV]'].values

# reference point comp_0
i0 = np.where(dirs == 'comp_0')[0][0]
E0, V0 = ener[i0], vols[i0]

dE = ener - E0
if delta_volume_units_percentage:
    dV = 100 * (vols - V0) / V0
    x_label = 'ΔV [%]'
else:
    dV = vols - V0
    x_label = 'ΔV [Å³]'

# read Sr-substituted data
sr_data = pd.read_csv('e-v-full_Sr.csv')
sr_dirs = sr_data['Directory'].values
sr_vols = sr_data['Cell volume [Å³]'].values
sr_ener = sr_data['Total energy [eV]'].values

i0_sr = np.where(sr_dirs == 'comp_0')[0][0]
Es0, Vs0 = sr_ener[i0_sr], sr_vols[i0_sr]

sr_dE = sr_ener - Es0
if delta_volume_units_percentage:
    sr_dV = 100 * (sr_vols - Vs0) / Vs0
else:
    sr_dV = sr_vols - Vs0

# plot setup
splot = SimplePlot()
splot.import_scatter_data_set(
    'Base', dV, dE, marker='o', color='C0'
)
splot.import_scatter_data_set(
    'Sr', sr_dV, sr_dE, marker='x', color='C1'
)

splot.ratio = (12, 6)
splot.x_label = x_label
splot.y_label = 'ΔE [eV]'
splot.show_grid = False
splot.show_legend = True

# explicitly add legend
plt.legend()

splot.store_plot = True
splot.store_path = './'
splot.show_plot = True

# draw
splot.plot(['Base', 'Sr'], plot_name='EOS.png')