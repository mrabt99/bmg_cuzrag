#!/usr/bin/env python3
"""Igual que step12_fig_compare_rates.py (mismos datos, paneles, colores y marcadores) pero con eje
vertical LOGARITMICO. Lee solo tabla_compare_rates_3tasas.csv; escribe figs/compare_rates_<panel>_log.png
(NO pisa las figuras lineales de step12).

Diferencias necesarias por el eje log:
  - Fraccion 0 (N=0) no existe en log: esos puntos se marcan con una "x" pegada a la base del eje y la linea
    se corta ahi (no se une un 0 con un valor positivo). Solo ocurre en Cu-BOO (15 de 51 puntos).
  - Barras de error Poisson sqrt(N): la barra inferior se detiene en la base del eje si pct - err <= base
    (cuando N es chico el error iguala o supera al valor).
  - Eje Y compartido dentro de cada grupo especie/definicion (como en step12), con margen de una fraccion de decada.
Marcador hueco = N<20 (igual que step12)."""
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

# limites Y por grupo: piso = decada/medio-decada por debajo del minimo positivo; techo = maximo (+err) con margen
ylim = {}
for sp, tc, _, _, g in PANELS:
    s = d[(d.especie == sp) & (d.tipo_cluster == tc)]
    pos = s[s.N > 0]
    lo, hi = pos.pct.min(), (pos.pct + pos.err).max()
    a, b = ylim.get(g, (np.inf, 0))
    ylim[g] = (min(a, lo), max(b, hi))
has_zero = {g: bool(((d.especie == sp) & (d.tipo_cluster == tc) & (d.N == 0)).any()) for sp, tc, _, _, g in PANELS}
# con ceros: una decada extra debajo del minimo positivo para que la banda de "x" quede separada de los puntos con N=1
ylim = {g: (10 ** (np.floor(np.log10(lo) * 2) / 2 - (0.6 if has_zero[g] else 0.15)), hi * 1.6) for g, (lo, hi) in ylim.items()}

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "text.color": INK, "figure.facecolor": SURF, "axes.facecolor": SURF})
for sp, tc, fname, title, g in PANELS:
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    s = d[(d.especie == sp) & (d.tipo_cluster == tc)]
    lo, hi = ylim[g]
    zero_y = lo * 1.3                                   # banda visual de los N=0 (justo sobre la base, separada de N>=1)
    any_zero = (s.N == 0).any()
    for ri, (rate, col, mk, lab) in enumerate(RATES):
        r = s[s.rate == rate].sort_values("T")
        assert len(r) == 17
        p = r[r.N > 0]
        lower = np.clip(p.pct - np.maximum(p.pct - p.err, lo), 0, None)      # la barra inferior no baja de la base
        ax.errorbar(p["T"], p.pct, yerr=[lower, p.err], color=col, lw=1.6, elinewidth=0.9, capsize=2, alpha=0.95, zorder=2)
        full, low = p[p.N >= 20], p[p.N < 20]
        ax.plot(full["T"], full.pct, ls="none", marker=mk, ms=5.5, mfc=col, mec=SURF, mew=0.8, color=col, zorder=3)
        ax.plot(low["T"], low.pct, ls="none", marker=mk, ms=5.5, mfc=SURF, mec=col, mew=1.3, color=col, zorder=3)
        z = r[r.N == 0]
        if len(z):
            # desplazo +-18 K solo para que las "x" de las 3 tasas no se tapen
            ax.plot(z["T"] + (ri - 1) * 18, np.full(len(z), zero_y), ls="none", marker="x", ms=6, mew=1.4, color=col, zorder=3)
    ax.set_yscale("log"); ax.set_ylim(lo, hi)
    ax.set_xlim(250, 1950); ax.set_xticks(range(300, 1901, 200))
    ax.set_xlabel("T (K)"); ax.set_ylabel(f"fracción de átomos de {sp} (%, escala log)")
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.grid(axis="y", which="major", color=GRID, lw=0.7); ax.set_axisbelow(True)
    for k in ("top", "right"): ax.spines[k].set_visible(False)
    hs = [Line2D([0], [0], color=c, lw=1.6, marker=m, ms=5.5, mfc=c, mec=SURF, label=l) for _, c, m, l in RATES]
    if (s.N < 20).any():
        hs.append(Line2D([0], [0], ls="none", marker="o", ms=5.5, mfc=SURF, mec=INK2, mew=1.3, label="N < 20"))
    if any_zero:
        hs.append(Line2D([0], [0], ls="none", marker="x", ms=6, mew=1.4, color=INK2, label="N = 0 (no representable en log)"))
    ax.legend(handles=hs, frameon=False, loc="upper right", fontsize=8.5, labelcolor=INK2)
    fig.text(0.012, 0.012, "Barras: error Poisson √N por punto de T (la inferior se corta en la base). Voronoi: OVITO default, índice exacto.",
             fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(f"figs/compare_rates_{fname}_log.png", dpi=150, facecolor=SURF)
    plt.close(fig)
    print(f"-> figs/compare_rates_{fname}_log.png  (Y {lo:.4g}-{hi:.4g} %, grupo {g})")
