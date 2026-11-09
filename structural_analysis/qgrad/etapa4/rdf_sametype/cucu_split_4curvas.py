#!/usr/bin/env python3
"""Desdoblamiento del 2do pico Cu-Cu para las 4 curvas. El script que genero CuCu_split_summary.csv original no esta en el
arbol; se RECUPERARON las ventanas por busqueda exhaustiva (unicas que reproducen exactamente r_peak1/r_dip/r_peak2 de las 51
filas originales): pico1 (4.0,4.9) A, valle (4.7,5.3) A, pico2 (5.0,5.7) A; depth = (g1+g2)/2 - g_dip (misma definicion que
ZrCu_split_summary). Se verifica que las 3 tasas originales se reproduzcan (assert) antes de escribir qgrad.
Salidas: CuCu_split_summary_4curvas.csv, ../fig_Q_CuCu_split_depth_4curvas.png"""
import csv
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
W1, WD, W2 = (4.0, 4.9), (4.7, 5.3), (5.0, 5.7)
g = {}
for r in csv.DictReader(open(HERE / "rdf_sametype_all.csv")):
    g.setdefault((r["rate"], int(float(r["T"]))), []).append((float(r["r"]), float(r["g_CuCu"])))
g = {k: np.array(v) for k, v in g.items()}
def ext(r, gg, w, f):
    m = (r > w[0]) & (r < w[1]); i = f(gg[m]); return r[m][i], gg[m][i]
rows = []
for (rate, T), a in sorted(g.items(), key=lambda kv: (kv[0][0], -kv[0][1])):
    r, gg = a[:, 0], a[:, 1]
    r1, g1 = ext(r, gg, W1, np.argmax); rd, gd = ext(r, gg, WD, np.argmin); r2, g2 = ext(r, gg, W2, np.argmax)
    rows.append([rate, T, r1, g1, rd, gd, r2, g2, (g1 + g2) / 2 - gd])
# verificacion contra el resumen original
orig = {(r["rate"], int(r["T"])): [float(r[c]) for c in ("r_peak1", "g_peak1", "r_dip", "g_dip", "r_peak2", "g_peak2", "depth")]
        for r in csv.DictReader(open(HERE.parent.parent.parent / "etapa4" / "rdf_sametype" / "CuCu_split_summary.csv"))}
mx = max(np.max(np.abs(np.array(r[2:]) - np.array(orig[(r[0], r[1])]))) for r in rows if r[0] != "qgrad")
print(f"reproduccion de CuCu_split_summary.csv original (51 filas): max|diff| = {mx:.2e}")
assert mx < 6e-5, "las ventanas no reproducen el resumen original"
with open(HERE / "CuCu_split_summary_4curvas.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "T", "r_peak1", "g_peak1", "r_dip", "g_dip", "r_peak2", "g_peak2", "depth"])
    for row in rows: w.writerow([row[0], row[1]] + [f"{v:.4f}" for v in row[2:]])
COL = {"1e12": "#d62728", "1e11": "#e6b800", "2e10": "#1f3fd6", "qgrad": "#7b3fbf"}
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
for rate in ["1e12", "1e11", "2e10", "qgrad"]:
    s = sorted([r for r in rows if r[0] == rate], key=lambda r: r[1])
    ax[0].plot([r[1] for r in s], [r[8] for r in s], "-o", ms=3.5, color=COL[rate], label=rate)
    ax[1].plot([r[1] for r in s], [r[6] for r in s], "-o", ms=3.5, color=COL[rate], label=rate)
ax[0].set(xlabel="T (K)", ylabel="profundidad del valle Cu-Cu", title="Profundidad del valle del 2do pico Cu-Cu")
ax[1].set(xlabel="T (K)", ylabel="r del 2do sub-pico (A)", title="Posicion del 2do sub-pico Cu-Cu")
for a in ax: a.grid(alpha=.25); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(HERE.parent / "fig_Q_CuCu_split_depth_4curvas.png", dpi=170)
for T in (1500, 900, 800, 700, 600, 500, 300):
    print(f"T={T:4d} depth: " + "  ".join(f"{k}={[r[8] for r in rows if r[0]==k and r[1]==T][0]:.3f}" for k in ["1e12", "1e11", "2e10", "qgrad"]) +
          "   r_pico2: " + "  ".join(f"{k}={[r[6] for r in rows if r[0]==k and r[1]==T][0]:.2f}" for k in ["2e10", "qgrad"]))
