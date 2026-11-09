#!/usr/bin/env python3
"""Curvas esfuerzo-deformacion de los barridos de la matriz (lee manifest.csv y los *.curve.dat de cada caso 'ok'; no modifica nada, solo escribe figuras nuevas en ../figs/).
Cada curva = media sobre las direcciones x, y, z del mismo vidrio (sigma del eje de traccion); banda = min-max entre direcciones (ruido termico + anisotropia del MISMO vidrio, no variacion entre vidrios).
Figuras (prefijo curvas_):
  curvas_tasa_def_300K.png   un panel por muestra (2e10, 1e11, 1e12, qgrad, qgrad_ext): colores = tasa de deformacion 1e8/1e9/1e10, 300 K
  curvas_T_1e9.png           un panel por muestra con barrido en T (2e10, 1e11, 1e12): colores = T (300-900 K), tasa 1e9
  curvas_historia_300K.png   un panel por tasa de deformacion; colores = muestra (historia termica); fila superior completa, inferior zoom al pico
  curvas_direcciones_300K_1e9.png   las 3 direcciones por separado para cada muestra (reproducibilidad), 300 K, 1e9
Uso (venv_structural): python stress_strain.py"""
import os, csv, glob
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); os.makedirs(f"{M}/figs", exist_ok=True)
CI = {"x": 2, "y": 3, "z": 4}
SAMPLES = ["2e10", "1e11", "1e12", "qgrad", "qgrad_ext"]
NAME = {"2e10": "2e10 K/s", "1e11": "1e11 K/s", "1e12": "1e12 K/s", "qgrad": "qgrad (graduado)", "qgrad_ext": "qgrad_ext (graduado, más lento)"}
COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad": "#7b4fc4", "qgrad_ext": "#b5338a"}
LS = {"2e10": "-", "1e11": "-", "1e12": "-", "qgrad": "--", "qgrad_ext": "--"}
INK, GRID = "#0b0b0b", "#e6e5e1"

curves = defaultdict(dict)   # (sample, T, erate) -> {dir: (strain, sigma)}
for r in csv.DictReader(open(f"{M}/manifest.csv")):
    if r["stage"] == "extra" or not r["status"].startswith("ok"): continue
    p = os.path.realpath(os.path.join(M, r["path"])); f = glob.glob(f"{p}/*.curve.dat")
    if not f: continue
    d = np.loadtxt(f[0], comments="#"); curves[(r["sample"], float(r["T"]), r["erate"])][r["dir"]] = (d[:, 1], d[:, CI[r["dir"]]])

def band(key):
    """media y min/max sobre direcciones; todas comparten la grilla de strain (curva promediada cada 0.001)"""
    c = curves[key]; n = min(len(v[0]) for v in c.values()); e = next(iter(c.values()))[0][:n]
    S = np.array([c[k][1][:n] for k in sorted(c)]); return e, S.mean(0), S.min(0), S.max(0), len(c)

def style(ax, xlim=(0, 0.2), ylim=None):
    ax.set(xlim=xlim); ax.grid(color=GRID, lw=.7); ax.set_axisbelow(True)
    if ylim: ax.set_ylim(*ylim)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

def draw(ax, key, color, label=None, ls="-", lw=1.5):
    e, m, lo, hi, n = band(key); ax.fill_between(e, lo, hi, color=color, alpha=.22, lw=0); ax.plot(e, m, color=color, lw=lw, ls=ls, label=label)

NOTE = "Línea: media sobre x, y, z del mismo vidrio; banda: mín–máx entre direcciones (no incluye variación entre vidrios)"
def seq(n, cmap="viridis"):
    cm = plt.get_cmap(cmap); return [cm(v) for v in np.linspace(0.05, 0.85, n)]

# 1) barrido de tasa de deformacion a 300 K, un panel por muestra
ER = ["1e8", "1e9", "1e10"]; ec = dict(zip(ER, seq(3, "plasma")))
fg, axs = plt.subplots(1, 5, figsize=(20, 4.4), sharey=True, constrained_layout=True)
for ax, s in zip(axs, SAMPLES):
    for er in ER:
        if (s, 300.0, er) in curves: draw(ax, (s, 300.0, er), ec[er], label=f"{er} s⁻¹")
    ax.set_title(NAME[s], color=COL[s], fontsize=11); ax.set_xlabel("deformación ε"); style(ax, ylim=(0, 3.1))
axs[0].set_ylabel("σ del eje de tracción (GPa)"); axs[0].legend(title="tasa de deformación", frameon=False, fontsize=9, loc="lower right")
fg.suptitle("Curvas esfuerzo–deformación: barrido de tasa de deformación (300 K)\n" + NOTE, fontsize=10.5); fg.savefig(f"{M}/figs/curvas_tasa_def_300K.png", dpi=140); plt.close(fg)

# 2) barrido de temperatura a 1e9, un panel por muestra (solo las que lo tienen)
TS = sorted({k[1] for k in curves if k[2] == "1e9"}); tc = dict(zip(TS, seq(len(TS), "coolwarm")))
tc = dict(zip(TS, seq(len(TS), "viridis")))
sm_T = [s for s in SAMPLES if sum(1 for k in curves if k[0] == s and k[2] == "1e9") > 1]
fg, axs = plt.subplots(1, len(sm_T), figsize=(5.6 * len(sm_T), 4.4), sharey=True, constrained_layout=True)
for ax, s in zip(np.atleast_1d(axs), sm_T):
    for T in TS:
        if (s, T, "1e9") in curves: draw(ax, (s, T, "1e9"), tc[T], label=f"{T:.0f} K")
    ax.set_title(NAME[s], color=COL[s], fontsize=11); ax.set_xlabel("deformación ε"); style(ax, ylim=(0, 3.1))
np.atleast_1d(axs)[0].set_ylabel("σ del eje de tracción (GPa)"); np.atleast_1d(axs)[0].legend(title="T", frameon=False, fontsize=9, loc="upper right")
fg.suptitle("Curvas esfuerzo–deformación: barrido de temperatura (tasa de deformación 1e9 s⁻¹; Tg ≈ 642, 666, 701 K para 2e10, 1e11, 1e12)\n" + NOTE, fontsize=10.5)
fg.savefig(f"{M}/figs/curvas_T_1e9.png", dpi=140); plt.close(fg)

# 3) efecto de la historia termica a 300 K, un panel por tasa de deformacion; fila de abajo = zoom al pico
fg, axs = plt.subplots(2, 3, figsize=(15, 8.2), constrained_layout=True)
for j, er in enumerate(ER):
    for i, (xl, yl) in enumerate([((0, 0.2), (0, 3.1)), ((0.03, 0.13), None)]):
        ax = axs[i, j]
        for s in SAMPLES:
            if (s, 300.0, er) in curves: draw(ax, (s, 300.0, er), COL[s], label=NAME[s] if (i == 0 and j == 2) else None, ls=LS[s], lw=1.4)
        ax.set_xlabel("deformación ε"); style(ax, xlim=xl)
        if i == 0: ax.set_title(f"tasa de deformación {er} s⁻¹", fontsize=11); ax.set_ylim(*yl)
        else:
            lo = max(band((s, 300.0, er))[1][(band((s, 300.0, er))[0] >= .03) & (band((s, 300.0, er))[0] <= .13)].min() for s in SAMPLES if (s, 300.0, er) in curves)
            hi = max(band((s, 300.0, er))[3][(band((s, 300.0, er))[0] <= .13)].max() for s in SAMPLES if (s, 300.0, er) in curves)
            ax.set_ylim(lo * 0.93, hi * 1.04); ax.set_title("zoom al pico", fontsize=10)
    axs[i, 0].set_ylabel("σ del eje de tracción (GPa)"); axs[1, 0].set_ylabel("σ del eje de tracción (GPa)")
axs[0, 2].legend(title="muestra (enfriamiento)", frameon=False, fontsize=9, loc="lower right")
fg.suptitle("Curvas esfuerzo–deformación: efecto de la historia térmica del vidrio (300 K)\n" + NOTE, fontsize=10.5)
fg.savefig(f"{M}/figs/curvas_historia_300K.png", dpi=140); plt.close(fg)

# 4) reproducibilidad: x, y, z por separado (300 K, 1e9), un panel por muestra
fg, axs = plt.subplots(1, 5, figsize=(20, 4.2), sharey=True, constrained_layout=True)
DC = {"x": "#2a78d6", "y": "#eb6834", "z": "#1baf7a"}
for ax, s in zip(axs, SAMPLES):
    for dname, (e, sg) in sorted(curves[(s, 300.0, "1e9")].items()): ax.plot(e, sg, color=DC[dname], lw=1.1, label=f"tracción en {dname}")
    ax.set_title(NAME[s], color=COL[s], fontsize=11); ax.set_xlabel("deformación ε"); style(ax, ylim=(0, 3.1))
axs[0].set_ylabel("σ del eje de tracción (GPa)"); axs[0].legend(frameon=False, fontsize=9, loc="lower right")
fg.suptitle("Curvas esfuerzo–deformación por dirección de tracción (300 K, 1e9 s⁻¹): mismo vidrio, tres direcciones", fontsize=10.5)
fg.savefig(f"{M}/figs/curvas_direcciones_300K_1e9.png", dpi=140); plt.close(fg)
print("curvas cargadas:", len(curves), "grupos;", sum(len(v) for v in curves.values()), "curvas; figuras en", f"{M}/figs/")
