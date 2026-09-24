#!/usr/bin/env python3
"""Tg_comparativo.png centrado en la zona de la transicion (T <= T_MAX).
Mismo calculo de Tg que run_tg_1e11.py (ventanas de TG_WINDOWS); solo
cambia el rango de T mostrado y la escala Y, que se ajusta a los datos
visibles. Ademas de los puntos binned se dibujan las dos rectas de ajuste
(rama liquida / rama vitrea) de cada tasa, extendidas hasta su interseccion.
No pisa Tg_comparativo.png: escribe Tg_comparativo_zoom.png.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compute_tg import load_thermo, bin_by_temperature, two_branch_tg, TG_WINDOWS

T_BIN = 10.0
T_MIN, T_MAX = 250.0, 1000.0
RATE_COLOR = {"2e10": "tab:blue", "1e11": "tab:orange", "1e12": "tab:red"}
RATE_CSV = {r: f"{r}/thermo_{r}.csv" for r in RATE_COLOR}

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ylims = {0: [np.inf, -np.inf], 1: [np.inf, -np.inf]}
for rate, color in RATE_COLOR.items():
    T, Ep, V = load_thermo(RATE_CSV[rate])
    Tb, Vb = bin_by_temperature(T, V, T_BIN)
    _, Epb = bin_by_temperature(T, Ep, T_BIN)
    sel = (Tb >= T_MIN) & (Tb <= T_MAX)
    for ax_i, (yb, key, lab) in enumerate([(Vb, "V", "Tg(V)"), (Epb, "Ep", "Tg(Ep)")]):
        Tg, (ml, cl), (mg, cg) = two_branch_tg(Tb, yb, TG_WINDOWS["liquid"], TG_WINDOWS["glass"])
        ax = axes[ax_i]
        ax.plot(Tb[sel], yb[sel], ".", ms=4, color=color, alpha=0.7)
        tl = np.array([TG_WINDOWS["liquid"][0], T_MAX])
        tg_ = np.array([T_MIN, TG_WINDOWS["glass"][1]])
        ax.plot(tl, ml * tl + cl, "-", lw=0.8, color=color, alpha=0.6)
        ax.plot(tg_, mg * tg_ + cg, "-", lw=0.8, color=color, alpha=0.6)
        ax.axvline(Tg, color=color, ls="--", lw=1, label=f"{rate}: {lab}={Tg:.0f}K")
        ax.plot(Tg, ml * Tg + cl, "x", color=color, ms=10, mew=2)
        ylims[ax_i][0] = min(ylims[ax_i][0], yb[sel].min())
        ylims[ax_i][1] = max(ylims[ax_i][1], yb[sel].max())

for ax_i, (title, ylabel) in enumerate([("Tg por volumen", "V (A^3)"),
                                        ("Tg por energia potencial", "Ep (eV/at)")]):
    ax = axes[ax_i]
    lo, hi = ylims[ax_i]
    pad = 0.04 * (hi - lo)
    ax.set_xlim(T_MIN, T_MAX)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_title(title); ax.set_xlabel("T (K)"); ax.set_ylabel(ylabel)
    ax.legend(loc="upper left")
fig.suptitle(f"Tg comparativo (zoom {T_MIN:.0f}-{T_MAX:.0f}K) -- ventanas: liquido 900-1200K, vitreo 300-600K")
fig.tight_layout()
fig.savefig("Tg_comparativo_zoom.png", dpi=130)
print("Tg_comparativo_zoom.png escrito")
