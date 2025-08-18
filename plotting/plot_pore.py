#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import numpy as np
import pandas as pd
from aim2dat.plots import SimplePlot

e_v = pd.read_csv('e-v-full.csv')
pore = pd.read_csv('pore_data.csv')
if 'Directory' not in e_v:
    raise KeyError("'Directory' missing in e-v-full.csv")

merged = pd.merge(
    e_v[['Directory', 'Cell volume [Å³]']],
    pore[[
        'Directory',
        'Largest_Included_Sphere',
        'Largest_Free_Sphere',
    ]],
    on='Directory',
    how='inner',
)

dirs = merged['Directory'].values
vols = merged['Cell volume [Å³]'].values
lis = merged['Largest_Included_Sphere'].values
lfs = merged['Largest_Free_Sphere'].values

idx0 = np.where(dirs == 'comp_0')[0]
if not idx0.size:
    raise ValueError("Reference 'comp_0' not found")
V0 = vols[idx0[0]]
dV = 100.0 * (vols - V0) / V0

plot = SimplePlot()
plot.backend = 'matplotlib'
plot.import_scatter_data_set('LIS', dV, lis, marker='o', color='teal')
plot.import_scatter_data_set('LFS', dV, lfs, marker='s', color='orange')
plot.ratio = (9, 6)
plot.x_label = 'ΔV [%]'
plot.y_label = 'Pore diameter [Å]'
plot.custom_xticks=np.arange(-4, 5, 1)
plot.x_range= (-4, 4)
plot.show_grid = False
plot.show_legend = False
plot.store_plot = False
plot.y_range = (6, 18)
plot.dpi         = 300


fig = plot.plot(['LIS', 'LFS'], plot_name=None)
ax = fig.axes[0]
names = ['LIS', 'LFS']
for line, name in zip(ax.lines, names):
    line.set_label(name)

ax.legend(loc='center right', fontsize=18)
fig.tight_layout()
fig.savefig('pore_data.png', dpi=300)
fig.show()
