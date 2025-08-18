#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator

# ==== USER CONFIGURATION ====
plot_all_minima = False
start_index = 1
delta_volume_units_percentage = True
# ============================

def read_helmholtz_volume(fname):
    data = {}
    temp = None
    vols, es = [], []
    with open(fname) as f:
        for line in f:
            if line.startswith("# Temperature:"):
                if temp is not None and temp % 100 == 0:
                    data[temp] = pd.DataFrame({
                        "Volume": vols,
                        "Free Energy": es
                    })
                temp = float(line.split(":", 1)[1])
                vols, es = [], []
            elif line.strip() and not line.startswith("#"):
                v, e = map(float, line.split())
                vols.append(v)
                es.append(e)
    if temp is not None and temp % 100 == 0:
        data[temp] = pd.DataFrame({
            "Volume": vols,
            "Free Energy": es
        })
    return data

def read_minima_points(fname):
    mins = []
    capture = False
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# Minimas"):
                capture = True
                continue
            if capture and line:
                try:
                    v, e = map(float, line.split())
                    mins.append((v, e))
                except ValueError:
                    pass
    return mins

def read_reference_volume(fname):
    vols = []
    with open(fname) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                v, _ = map(float, line.split())
                vols.append(v)
    if not vols:
        raise ValueError("No volumes in e-v.dat")
    return vols[len(vols) // 2]

def adjust_free_energy_reference(data, ref_temp=300):
    if ref_temp not in data:
        raise ValueError(f"{ref_temp} K not in dataset")
    e0 = data[ref_temp]["Free Energy"].min()
    for df in data.values():
        df["Free Energy"] -= e0
    return data

def shift_volumes(data, ref_vol, percentage=False):
    for df in data.values():
        if percentage:
            df["Volume"] = 100 * (df["Volume"] - ref_vol) / ref_vol
        else:
            df["Volume"] -= ref_vol
    return data

def shift_minima_volumes(mins, ref_vol, percentage=False):
    if percentage:
        return [(100*(v-ref_vol)/ref_vol, e) for v, e in mins]
    return [(v-ref_vol, e) for v, e in mins]

def plot_helmholtz_volume(data, minima, percentage=False, fname="helmholtz_plot.png"):
    fig, ax = plt.subplots(figsize=(6, 9))
    temps = sorted(data)
    colors = plt.cm.plasma(np.linspace(0, 1, len(temps)))
    for temp, c in zip(temps, colors):
        df = data[temp]
        ax.plot(df["Volume"], df["Free Energy"], "o-",
                label=f"T={int(temp)} K", color=c)
    if minima:
        vs, es = zip(*minima)
        ax.plot(vs, es, "-", lw=2, color="red", zorder=1)
        ax.scatter(vs, es, color="red", s=50, label="Minima", zorder=2)

    # auto‐set limits and ticks
    all_vol = np.hstack([df["Volume"].values for df in data.values()])
    all_fe  = np.hstack([df["Free Energy"].values for df in data.values()])
    x_min, x_max = all_vol.min(), all_vol.max()
    y_min, y_max = all_fe.min(), all_fe.max()-1
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xticks(AutoLocator().tick_values(x_min, x_max))
    ax.set_yticks(AutoLocator().tick_values(y_min, y_max))
    ax.tick_params(axis='both', labelsize=17)


    xlabel = r"$\Delta V$ [%]" if percentage else r"$\Delta V$ / Å³"
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel("Free Energy / eV", fontsize=19)
    ax.legend(loc="lower right", fontsize=15)
    fig.tight_layout()
    fig.savefig(fname, dpi=300)
    plt.show()

def main():
    data = read_helmholtz_volume("helmholtz-volume.dat")
    data = adjust_free_energy_reference(data, ref_temp=300)
    all_min = read_minima_points("helmholtz-volume_fitted.dat")
    ref_vol = read_reference_volume("e-v.dat")
    data = shift_volumes(data, ref_vol,
                         percentage=delta_volume_units_percentage)
    all_min = shift_minima_volumes(all_min, ref_vol,
                                   percentage=delta_volume_units_percentage)
    minima = all_min if plot_all_minima else all_min[start_index - 1:]
    plot_helmholtz_volume(data, minima,
                          percentage=delta_volume_units_percentage)

if __name__ == "__main__":
    main()
