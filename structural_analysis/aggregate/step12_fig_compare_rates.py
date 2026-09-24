#!/usr/bin/env python3
"""Tarea R (paso 2) -- compare_rates: fraccion (%) vs T, 3 tasas superpuestas, un PNG por
tipo de cluster, SIN agregar por regimen. Lee solo tabla_compare_rates_3tasas.csv.
Barras = error Poisson sqrt(N) por punto de T (en % de la especie). Marcador hueco = N<20.
Eje Y compartido dentro de cada grupo especie/definicion; eje T comun en todos los paneles.
Paleta categorica en orden fijo (blue/orange/aqua), validada (CVD ok; aqua <3:1 -> etiquetas + tabla)."""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

NSP = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
RATES = [("2e10", "#2a78d6", "o", "2e10 K/s"), ("1e11", "#eb6834", "s", "1e11 K/s"), ("1e12", "#1baf7a", "^", "1e12 K/s")]
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
PANELS = [  # (especie, tipo_cluster, archivo, titulo, grupo-Y)
    ("Cu", "ico_boo",       "Cu_ico_BOO",      "Cu — icosaedro (BOO: q6 ≥ 0.6, Ŵ6 ≤ −0.16)", "Cu_boo"),
    ("Cu", "ico_voronoi",   "Cu_ico_Voronoi",  "Cu — icosaedro Voronoi ⟨0,0,12,0⟩", "Cu_vor"),
    ("Zr", "Z14_voronoi",   "Zr_Z14",          "Zr — Z14, Voronoi ⟨0,0,12,2⟩", "Zr_fk"),
    ("Zr", "Z15_voronoi",   "Zr_Z15",          "Zr — Z15, Voronoi ⟨0,0,12,3⟩", "Zr_fk"),
    ("Zr", "Z16_voronoi",   "Zr_Z16",          "Zr — Z16, Voronoi ⟨0,0,12,4⟩ (límite superior)", "Zr_fk"),
    ("Ag", "ico_voronoi",   "Ag_ico_Voronoi",  "Ag — icosaedro Voronoi ⟨0,0,12,0⟩", "Ag_vor"),
]
d = pd.read_csv("tabla_compare_rates_3tasas.csv", dtype={"rate": str})
d["pct"] = 100 * d.fraccion
d["err"] = 100 * np.sqrt(d.N) / d.especie.map(NSP)

ymax = {}
for sp, tc, _, _, g in PANELS:
    s = d[(d.especie == sp) & (d.tipo_cluster == tc)]
    ymax[g] = max(ymax.get(g, 0), (s.pct + s.err).max())
ymax = {g: v * 1.08 for g, v in ymax.items()}
# techo Y propio de un panel (rompe el eje compartido del grupo a proposito): Z16 (max 1.55 % con barra) no necesita el
# techo de Z15 (2.8 %); 1.75 % lo pidio el usuario. Ver auditoria: no recorta datos ni barras.
YMAX_OVERRIDE = {"Zr_Z16": 1.75}
# uso: python step12_fig_compare_rates.py [archivo ...]  -> regenera solo esos paneles (p.ej. Zr_Z16); sin args, todos
ONLY = set(sys.argv[1:])

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "text.color": INK, "figure.facecolor": SURF, "axes.facecolor": SURF})
for sp, tc, fname, title, g in PANELS:
    if ONLY and fname not in ONLY:
        continue
    ytop = YMAX_OVERRIDE.get(fname, ymax[g])
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    s = d[(d.especie == sp) & (d.tipo_cluster == tc)]
    for rate, col, mk, lab in RATES:
        r = s[s.rate == rate].sort_values("T")
        assert len(r) == 17
        ax.errorbar(r["T"], r.pct, yerr=r.err, color=col, lw=1.6, elinewidth=0.9, capsize=2, alpha=0.95, zorder=2)
        full, low = r[r.N >= 20], r[r.N < 20]
        ax.plot(full["T"], full.pct, ls="none", marker=mk, ms=5.5, mfc=col, mec=SURF, mew=0.8, color=col, zorder=3)
        ax.plot(low["T"], low.pct, ls="none", marker=mk, ms=5.5, mfc=SURF, mec=col, mew=1.3, color=col, zorder=3)
    ax.set_xlim(250, 1950); ax.set_xticks(range(300, 1901, 200)); ax.set_ylim(0, ytop)
    ax.set_xlabel("T (K)"); ax.set_ylabel(f"fracción de átomos de {sp} (%)")
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.7); ax.set_axisbelow(True)
    for k in ("top", "right"): ax.spines[k].set_visible(False)
    hs = [Line2D([0], [0], color=c, lw=1.6, marker=m, ms=5.5, mfc=c, mec=SURF, label=l) for _, c, m, l in RATES]
    if (s.N < 20).any():
        hs.append(Line2D([0], [0], ls="none", marker="o", ms=5.5, mfc=SURF, mec=INK2, mew=1.3, label="N < 20"))
    ax.legend(handles=hs, frameon=False, loc="upper right", fontsize=8.5, labelcolor=INK2)
    fig.text(0.012, 0.012, "Barras: error Poisson √N por punto de T. Voronoi: OVITO default, índice exacto.", fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(f"figs/compare_rates_{fname}.png", dpi=150, facecolor=SURF)
    plt.close(fig)
    print(f"-> figs/compare_rates_{fname}.png  (Y 0-{ytop:.2f} %, grupo {g})")
