#!/usr/bin/env python3
"""Tendencias con la 1e12 UNIFICADA (1e12u, etapas 6 y 7; 17/11/2026). Script nuevo: no toca trends.py ni tendencias.csv.
Lee ../resumen.csv. Agrega por (muestra, T, erate) la media y el desvio (ddof=1) sobre x, y, z de E, sigma_max, eps(sigma_max),
sigma_flujo, delta_sigma, psi(0.02-0.04) y psi(0.04-0.06) -> ../tendencias_1e12u.csv (muestras 2e10, 1e11, 1e12u, 1e12 vieja, qgrad_ext).
Imprime: 1e12u contra 1e12 vieja caso por caso (mismo T y tasa), la brecha 2e10 - 1e12u, m = d ln sigma_max / d ln erate a 300 K
y d sigma_max / dT entre 300 y 540 K. qgrad queda afuera (pedido del usuario).
Barras = desvio entre direcciones del MISMO vidrio (no incluye variacion entre vidrios independientes).
Uso: /scratch/negishi/aaabt/venv_structural/bin/python trends_1e12u.py"""
import os, csv, math
from collections import defaultdict
import numpy as np
M = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
Q = ["E_GPa", "sigma_max_GPa", "eps_sigma_max", "sigma_flujo_GPa", "delta_sigma_GPa", "psi_w02_04", "psi_w04_06"]
SAMPLES = ["2e10", "1e11", "1e12u", "1e12", "qgrad_ext"]
g = defaultdict(list)
for r in csv.DictReader(open(f"{M}/resumen.csv")):
    if r["sample"] in SAMPLES: g[(r["sample"], float(r["T"]), r["erate"])].append(r)
agg = {}
for k, rs in sorted(g.items(), key=lambda kv: (SAMPLES.index(kv[0][0]), kv[0][1], float(kv[0][2]))):
    agg[k] = {q: (np.mean([float(r[q]) for r in rs]), np.std([float(r[q]) for r in rs], ddof=1) if len(rs) > 1 else float("nan")) for q in Q}
    agg[k]["n"] = len(rs)
with open(f"{M}/tendencias_1e12u.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample", "T", "erate", "n"] + [f"{q}_{s}" for q in Q for s in ("media", "sd")])
    for (s, T, er), a in agg.items(): w.writerow([s, f"{T:g}", er, a["n"]] + [f"{a[q][i]:.5g}" for q in Q for i in (0, 1)])

print("1e12u vs 1e12 vieja (media ± sd x,y,z):")
for (s, T, er), a in agg.items():
    if s != "1e12u" or ("1e12", T, er) not in agg: continue
    b = agg[("1e12", T, er)]
    print(f"  {T:4.0f} K {er:5s} | σmax {a['sigma_max_GPa'][0]:.3f}±{a['sigma_max_GPa'][1]:.3f} (vieja {b['sigma_max_GPa'][0]:.3f}±{b['sigma_max_GPa'][1]:.3f})"
          f" | Δσ {a['delta_sigma_GPa'][0]:.3f} ({b['delta_sigma_GPa'][0]:.3f}) | E {a['E_GPa'][0]:.1f} ({b['E_GPa'][0]:.1f})"
          f" | σflujo {a['sigma_flujo_GPa'][0]:.3f} ({b['sigma_flujo_GPa'][0]:.3f}) | ε(σmax) {a['eps_sigma_max'][0]:.4f} ({b['eps_sigma_max'][0]:.4f})"
          f" | ψ0.02-0.04 {a['psi_w02_04'][0]:.4f} ({b['psi_w02_04'][0]:.4f})")
print("Brecha 2e10 - 1e12u en σmax y Δσ (GPa; % de 2e10):")
for (s, T, er), a in agg.items():
    if s != "2e10" or ("1e12u", T, er) not in agg: continue
    b = agg[("1e12u", T, er)]; d = a["sigma_max_GPa"][0] - b["sigma_max_GPa"][0]
    print(f"  {T:4.0f} K {er:5s}: σmax {d:.3f} ({100 * d / a['sigma_max_GPa'][0]:.1f} %), Δσ {a['delta_sigma_GPa'][0] - b['delta_sigma_GPa'][0]:.3f},"
          f" Δσ(1e12u)/Δσ(2e10) = {b['delta_sigma_GPa'][0] / a['delta_sigma_GPa'][0]:.2f}, σflujo {100 * (b['sigma_flujo_GPa'][0] / a['sigma_flujo_GPa'][0] - 1):+.1f} %")
print("m = d ln σmax / d ln erate a 300 K (1e8-1e9, 1e9-1e10, ajuste de las 3):")
for s in SAMPLES:
    try: y = [agg[(s, 300.0, e)]["sigma_max_GPa"][0] for e in ("1e8", "1e9", "1e10")]
    except KeyError: continue
    lx = np.log([1e8, 1e9, 1e10]); ly = np.log(y)
    print(f"  {s:9s}: {(ly[1] - ly[0]) / (lx[1] - lx[0]):.3f}, {(ly[2] - ly[1]) / (lx[2] - lx[1]):.3f}, {np.polyfit(lx, ly, 1)[0]:.3f}")
print("dσmax/dT y dE/dT entre 300 y 540 K (1e9; ajuste de 300, 450, 540):")
for s in SAMPLES:
    try: pts = [(T, agg[(s, T, "1e9")]) for T in (300.0, 450.0, 540.0)]
    except KeyError: continue
    Ts = [p[0] for p in pts]
    print(f"  {s:9s}: σmax {1e3 * np.polyfit(Ts, [p[1]['sigma_max_GPa'][0] for p in pts], 1)[0]:.2f} MPa/K, E {1e3 * np.polyfit(Ts, [p[1]['E_GPa'][0] for p in pts], 1)[0]:.1f} MPa/K")
print("->", f"{M}/tendencias_1e12u.csv")
