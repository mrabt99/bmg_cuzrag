#!/usr/bin/env python3
"""Tendencias de la matriz (lee ../resumen.csv). Agrega por (muestra, T, erate): media y desvio (ddof=1) sobre las direcciones x, y, z de E, sigma_max, eps(sigma_max), sigma_flujo,
delta_sigma, psi(0.02-0.04), psi(0.04-0.06), psi_max -> ../tendencias.csv. Figuras: ../figs/tendencias_T.png (a 1e9, vs T) y ../figs/tendencias_erate.png (a 300 K, vs tasa de deformacion).
Barras = desvio entre direcciones (MISMO vidrio): NO incluyen la variacion entre vidrios independientes.
Tambien imprime la sensibilidad a la tasa m = d ln(sigma_max)/d ln(erate) a 300 K y d(sigma_max)/dT entre 300 y 540 K por muestra.
Uso: python trends.py"""
import os, csv, math
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); os.makedirs(f"{M}/figs", exist_ok=True)
Q = [("E_GPa", "E (GPa)"), ("sigma_max_GPa", "σmax (GPa)"), ("eps_sigma_max", "ε(σmax)"), ("sigma_flujo_GPa", "σ flujo (GPa)"), ("delta_sigma_GPa", "Δσ = σmax − σflujo (GPa)"),
     ("psi_w02_04", "ψ(0.02–0.04)"), ("psi_w04_06", "ψ(0.04–0.06)"), ("psi_max", "ψ máx")]
g = defaultdict(list)
for r in csv.DictReader(open(f"{M}/resumen.csv")): g[(r["sample"], float(r["T"]), r["erate"])].append(r)
def ms(rs, c):
    v = np.array([float(r[c]) for r in rs]); return v.mean(), (v.std(ddof=1) if len(v) > 1 else float("nan"))
with open(f"{M}/tendencias.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample", "T", "erate", "n_dirs"] + [x for c, _ in Q for x in (c + "_media", c + "_sd")])
    for k in sorted(g, key=lambda k: (k[0], k[2], k[1])):
        row = [k[0], int(k[1]), k[2], len(g[k])]
        for c, _ in Q: m, s = ms(g[k], c); row += [f"{m:.5g}", f"{s:.3g}"]
        w.writerow(row)
COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a"}; MK = {"2e10": "o", "1e11": "s", "1e12": "^"}; SAMPLES = ["2e10", "1e11", "1e12"]
PAN = [("E_GPa", "E (GPa)"), ("sigma_max_GPa", "σmax (GPa)"), ("sigma_flujo_GPa", "σ flujo, ε 0.15–0.20 (GPa)"), ("delta_sigma_GPa", "Δσ (GPa)"), ("psi_w02_04", "ψ(0.02–0.04)"), ("eps_sigma_max", "ε(σmax)")]
def fig(xkey, xs_of, xlabel, fname, title, logx=False):
    fg, axs = plt.subplots(2, 3, figsize=(14, 7.2), constrained_layout=True)
    for ax, (c, nm) in zip(axs.flat, PAN):
        for s in SAMPLES:
            pts = [(xs_of(k), *ms(g[k], c)) for k in g if k[0] == s and xkey(k)]
            if not pts: continue
            pts.sort(); x, m, sd = zip(*pts)
            ax.errorbar(x, m, yerr=np.nan_to_num(sd), color=COL[s], marker=MK[s], ms=5.5, lw=1.5, capsize=3, label=s)
        ax.set(xlabel=xlabel, ylabel=nm); ax.grid(color="#e6e5e1", lw=.7); ax.set_axisbelow(True)
        if logx: ax.set_xscale("log")
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    axs.flat[0].legend(title="enfriamiento (K/s)", frameon=False, fontsize=9)
    fg.suptitle(title + "\nBarras: desvío entre x, y, z del mismo vidrio (no incluye variación entre vidrios)", fontsize=10.5); fg.savefig(f"{M}/figs/{fname}", dpi=140); plt.close(fg)
fig(lambda k: k[2] == "1e9", lambda k: k[1], "T (K)", "tendencias_T.png", "Dependencia con la temperatura (tasa de deformación 1e9 s⁻¹)")
fig(lambda k: k[1] == 300.0, lambda k: float(k[2]), "tasa de deformación (s⁻¹)", "tendencias_erate.png", "Dependencia con la tasa de deformación (300 K)", logx=True)
print("== sensibilidad a la tasa a 300 K: m = d ln(sigma_max) / d ln(erate) (ajuste con las tasas disponibles)")
for s in SAMPLES:
    pts = sorted((float(k[2]), ms(g[k], "sigma_max_GPa")[0]) for k in g if k[0] == s and k[1] == 300.0)
    if len(pts) >= 2: x = np.log([p[0] for p in pts]); y = np.log([p[1] for p in pts]); print(f"  {s}: m = {np.polyfit(x, y, 1)[0]:.4f}  (tasas {[p[0] for p in pts]}, sigma_max {[round(p[1],3) for p in pts]})")
print("== sigma_max y E vs T a 1e9 (media sobre x,y,z); pendiente lineal 300-540 K")
for s in SAMPLES:
    pts = sorted((k[1], ms(g[k], "sigma_max_GPa")[0], ms(g[k], "E_GPa")[0]) for k in g if k[0] == s and k[2] == "1e9")
    if len(pts) >= 2:
        sel = [p for p in pts if p[0] <= 540]; T = [p[0] for p in sel]
        print(f"  {s}: " + "  ".join(f"T={p[0]:.0f}: smax={p[1]:.3f} E={p[2]:.1f}" for p in pts) + (f"  | d(smax)/dT = {np.polyfit(T,[p[1] for p in sel],1)[0]*1000:.3f} MPa/K, dE/dT = {np.polyfit(T,[p[2] for p in sel],1)[0]*1000:.1f} MPa/K" if len(sel) >= 2 else ""))
print(f"tendencias.csv: {len(g)} grupos; figuras en {M}/figs/")
