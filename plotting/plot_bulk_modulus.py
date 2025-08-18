#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import numpy as np
import matplotlib.pyplot as plt

def plot_bulk_modulus(file_path, output_file="bulk_modulus_plot.png"):
    data = np.loadtxt(file_path)
    temperature = data[:, 0]
    bulk_modulus = data[:, 1]

    bulk_modulus_0 = bulk_modulus[0]

    plt.figure(figsize=(4, 6))
    plt.plot(temperature, bulk_modulus, marker='o', linestyle='-', color='b', label='Bulk Modulus')

    plt.xlabel("Temperature (K)")
    plt.ylabel("Bulk Modulus (GPa)")
    plt.legend()

    ax = plt.gca()
    yticks = list(ax.get_yticks())

    if bulk_modulus_0 not in yticks:
        yticks.append(bulk_modulus_0)
        yticks = sorted(yticks)

    ax.set_yticks(yticks)

    ytick_labels = [f"{ytick:.1f} " if abs(ytick - bulk_modulus_0) < 1e-3 else f"{ytick:.1f}" for ytick in yticks]
    ax.set_yticklabels(ytick_labels)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)

plot_bulk_modulus("bulk_modulus-temperature.dat")
