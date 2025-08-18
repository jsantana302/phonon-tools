#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import numpy as np
import pandas as pd
from matplotlib.ticker import AutoLocator
from aim2dat.plots import SimplePlot

delta_volume_units_percentage = True

# load data
e_v_data = pd.read_csv('e-v-full.csv')
dirs    = e_v_data['Directory'].values
vols    = e_v_data['Cell volume [Å³]'].values
energies= e_v_data['Total energy [eV]'].values

# reference point
i0 = np.where(dirs == 'comp_0')[0][0]
E0, V0 = energies[i0], vols[i0]

dE = energies - E0
dV = vols    - V0

if delta_volume_units_percentage:
    dV = 100 * dV / V0
    x_label = 'ΔV [%]'
else:
    x_label = 'ΔV [Å³]'

# compute auto limits & ticks
x_min, x_max = dV.min(), dV.max()
y_min, y_max = dE.min(), dE.max()

locator = AutoLocator()
x_ticks = locator.tick_values(x_min, x_max)
y_ticks = locator.tick_values(y_min, y_max)

# set up plot
splot = SimplePlot()
splot.backend       = 'matplotlib'
splot.import_scatter_data_set(
    'ΔE vs ΔV', dV, dE,
    marker='o', color='C0'
)

splot.ratio         = (9, 6)
splot.x_label       = x_label
splot.y_label       = 'ΔE [eV]'
splot.show_grid     = False
splot.show_legend   = True

# apply auto limits & ticks
splot.x_range       = (x_min, x_max)
splot.custom_xticks = x_ticks
splot.y_range       = (y_min, y_max)
splot.custom_yticks = y_ticks

splot.store_plot    = True
splot.store_path    = './'
splot.show_plot     = True
splot.dpi           = 300

# draw
splot.plot(
    ['ΔE vs ΔV'],
    plot_name='EOS.png'
)
