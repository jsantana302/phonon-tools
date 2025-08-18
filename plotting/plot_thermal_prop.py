#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from aim2dat.plots import SimplePlot

PANEL_A = [
    ('thermal_expansion.dat',
     'Thermal Expansion (×10⁻⁶ K⁻¹)',
     'thermal.png', 'C0'),
    ('Cp-temperature_polyfit.dat',
     'Heat Capacity (J/mol·K)',
     'cp.png', 'C1'),
    ('bulk_modulus-temperature.dat',
     'Bulk Modulus (GPa)',
     'bulk.png', 'C2'),
    ('gruneisen-temperature.dat',
     'Grüneisen Parameter',
     'gruneisen.png', 'C3'),
]
PANEL_A_OUT = 'panel_thermodynamic.png'

PANEL_B = [
    ('volume-temperature.dat',
     'Volume (×10² Å³)',
     'volume.png', 'C0'),
    ('thermal_expansion.dat',
     'Thermal Expansion (×10⁻⁶ K⁻¹)',
     'thermal.png', 'C1'),
    ('Cp-temperature_polyfit.dat',
     'Heat Capacity (J/mol·K)',
     'cp.png', 'C2'),
    ('bulk_modulus-temperature.dat',
     'Bulk Modulus (GPa)',
     'bulk.png', 'C3'),
    ('gruneisen-temperature.dat',
     'Grüneisen Parameter',
     'gruneisen.png', 'C4'),
    ('gibbs-temperature.dat',
     'Gibbs Free Energy (eV)',
     'gibbs.png', 'C5'),
]
PANEL_B_OUT = 'panel_all.png'

LABELS = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

def make_plot(data_file, y_label, out_png, color):
    T, Y = np.loadtxt(data_file, unpack=True)
    if data_file == 'volume-temperature.dat':
        Y /= 1e2
        y_label = 'Volume (×10² Å³)'
    if data_file == 'thermal_expansion.dat':
        Y *= 1e6
        y_label = 'Thermal Expansion (×10⁻⁶ K⁻¹)'
    sp = SimplePlot()
    sp.import_scatter_data_set('data', T, Y, marker='o', color=color)
    sp.ratio = (8, 5)
    sp.dpi = 300
    sp.subplot_nrows = 1
    sp.subplot_ncols = 1
    sp.subplot_tight_layout = True
    sp.x_label = 'Temperature (K)'
    sp.y_label = y_label
    sp.store_plot = True
    sp.store_path = './'
    sp.show_plot = False
    sp.plot(['data'], plot_title='', plot_name=out_png)

def assemble_panel(panel, out_png, nrows, ncols):
    # more horizontal: width scale=8, height scale=5
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(8*ncols, 5*nrows))
    axes = axes.flatten()
    for ax, (_, _, png, _), lbl in zip(axes, panel, LABELS):
        ax.imshow(mpimg.imread(png))
        ax.axis('off')
        ax.text(-0.04, 1.02, lbl,
                transform=ax.transAxes, ha='left',
                va='top', fontsize=19)
    plt.subplots_adjust(
        left=0.05, right=0.95,
        top=0.95, bottom=0.05,
        wspace=0.1, hspace=0.1
    )
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    for _, _, png, _ in panel:
        if os.path.exists(png):
            os.remove(png)
    print(f'✅ Panel saved as {out_png}')

def main():
    for f, y, out, c in PANEL_A:
        make_plot(f, y, out, c)
    assemble_panel(PANEL_A, PANEL_A_OUT, 2, 2)
    for f, y, out, c in PANEL_B:
        make_plot(f, y, out, c)
    assemble_panel(PANEL_B, PANEL_B_OUT, 2, 3)

if __name__ == '__main__':
    main()
