"""Zr, 2e10 K/s: fraccion (%) vs T de Z14/Z15/Z16 segun tres definiciones.
Evidencia del informe: la figura compare_n_Z1x.png (CN==k + umbral de cara)
no mide lo mismo que el indice exacto <0,0,12,k>. Eje Y log, mismo rango en los 3 paneles."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
f = pd.read_csv("a2_recalculo_51configs.csv", dtype={"rate": str})
a = pd.read_csv("../abundancia_2e10.csv")
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
C = {"ex": "#2a78d6", "cnd": "#eb6834", "cnt": "#1baf7a"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "text.color": INK, "figure.facecolor": SURF, "axes.facecolor": SURF})
fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), sharey=True)
for ax, z in zip(axes, ["Z14", "Z15", "Z16"]):
    s = f[(f.rate == "2e10") & (f.especie == "Zr") & (f.tc == z)].sort_values("T")
    t = a[a.especie == "Zr"].sort_values("T")
    ax.plot(s["T"], 100 * s.exact_cols25, "-o", ms=4, lw=1.6, color=C["ex"], mec=SURF, label="índice exacto ⟨0,0,12,k⟩ (OVITO default)")
    ax.plot(s["T"], 100 * s.cn_default, "-s", ms=4, lw=1.6, color=C["cnd"], mec=SURF, label="CN == 12+k, OVITO default")
    ax.plot(t["T"], 100 * t[f"n_{z}"], "-^", ms=4.5, lw=1.6, color=C["cnt"], mec=SURF, label="CN == 12+k + umbral de cara (figura compare_n_Z*)")
    ax.set_yscale("log"); ax.set_ylim(0.03, 60); ax.invert_xaxis()
    ax.set_title(f"Zr — {z}", loc="left", fontsize=10); ax.set_xlabel("T (K)")
    ax.grid(axis="y", color=GRID, lw=0.7); ax.set_axisbelow(True)
    for k in ("top", "right"): ax.spines[k].set_visible(False)
axes[0].set_ylabel("% de átomos de Zr (log)")
axes[0].legend(frameon=False, fontsize=7.6, loc="lower left", labelcolor=INK2)
fig.text(0.012, 0.012, "2e10 K/s, 17 configuraciones. Solo Zr: en Cu y Ag el índice exacto Z14/15/16 tiene N<20 casi en todo T.", fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig("fig_zr_definiciones_2e10.png", dpi=150, facecolor=SURF)
