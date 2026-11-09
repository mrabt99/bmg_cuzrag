#!/usr/bin/env python3
"""Fraccion (%) vs T de 5 curvas (2e10, 1e11, 1e12, qgrad, qgrad_ext) para Cu-ico (BOO y Voronoi) y Zr Z14/Z15/Z16, T 300-1000 K.
Barras: Poisson sqrt(N). Marcador hueco: N<20. 300 K de qgrad y qgrad_ext = 17r. Lee tabla_compare_rates_5curvas.csv. Salida figs/compare_ext_5curvas.png"""
import csv, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
NSP = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
R = [("2e10", "#2a78d6", "o", "2e10"), ("1e11", "#eb6834", "s", "1e11"), ("1e12", "#1baf7a", "^", "1e12"), ("qgrad", "#7b3fbf", "D", "qgrad"), ("qgrad_ext", "#0b0b0b", "P", "qgrad_ext (lenta hasta 500 K)")]
SURF, INK2, GRID = "#fcfcfb", "#52514e", "#e6e5e1"
PAN = [("Cu", "ico_boo", "Cu-ico (BOO)"), ("Cu", "ico_voronoi", "Cu-ico (Voronoi)"), ("Zr", "Z14_voronoi", "Zr Z14"), ("Zr", "Z15_voronoi", "Zr Z15"), ("Zr", "Z16_voronoi", "Zr Z16"), ("Ag", "ico_voronoi", "Ag-ico (Voronoi)")]
D = {}
for r in csv.DictReader(open("tabla_compare_rates_5curvas.csv")):
    T = int(r["T"])
    if T == 300 and r["rate"] in ("qgrad", "qgrad_ext") and r["config"] != "17r": continue
    D[(r["rate"], r["especie"], r["tipo_cluster"], T)] = int(r["N"])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF})
fig, axs = plt.subplots(2, 3, figsize=(13, 7.2))
for ax, (sp, tc, ttl) in zip(axs.flat, PAN):
    for rate, col, mk, lab in R:
        Ts = sorted({k[3] for k in D if k[0] == rate and k[1] == sp and k[2] == tc and k[3] <= 1000})
        n = np.array([D[(rate, sp, tc, T)] for T in Ts]); y = 100 * n / NSP[sp]; e = 100 * np.sqrt(n) / NSP[sp]
        lw = 2.0 if rate == "qgrad_ext" else 1.3
        ax.errorbar(Ts, y, yerr=e, color=col, lw=lw, elinewidth=.8, capsize=1.5, alpha=.95, zorder=4 if rate == "qgrad_ext" else 2)
        ok = n >= 20; Ts = np.array(Ts)
        ax.plot(Ts[ok], y[ok], ls="none", marker=mk, ms=5, mfc=col, mec=SURF, mew=.7, zorder=5)
        ax.plot(Ts[~ok], y[~ok], ls="none", marker=mk, ms=5, mfc=SURF, mec=col, mew=1.2, zorder=5)
    ax.set_xlim(250, 1000); ax.set_ylim(bottom=0); ax.set_title(ttl, loc="left", fontsize=10); ax.set_xlabel("T (K)"); ax.set_ylabel(f"% de átomos de {sp}")
    ax.grid(axis="y", color=GRID, lw=.7); ax.set_axisbelow(True)
    for k in ("top", "right"): ax.spines[k].set_visible(False)
h, l = axs.flat[0].get_legend_handles_labels()
from matplotlib.lines import Line2D
fig.legend([Line2D([0], [0], color=c, lw=1.6, marker=m, ms=5, mfc=c, mec=SURF) for _, c, m, _ in R], [x[3] for x in R], loc="lower center", ncol=5, frameon=False, fontsize=9, labelcolor=INK2)
fig.text(.01, .005, "Barras: Poisson √N. Hueco: N<20. 300 K de qgrad y qgrad_ext = 17r (relajado). Voronoi OVITO default.", fontsize=7, color=INK2)
fig.tight_layout(rect=(0, .05, 1, 1)); fig.savefig("figs/compare_ext_5curvas.png", dpi=150, facecolor=SURF); print("-> figs/compare_ext_5curvas.png")
