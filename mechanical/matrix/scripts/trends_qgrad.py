#!/usr/bin/env python3
"""Tendencias a 300 K incluyendo qgrad y qgrad_ext (lee ../resumen.csv; NO modifica trends.py ni sus salidas).
Escribe ../tendencias_con_qgrad.csv (media y desvio sobre x,y,z por muestra/erate a 300 K) y ../figs/tendencias_erate_con_qgrad.png (mismos 6 paneles que tendencias_erate.png con 5 muestras).
Imprime: sigma_max a 1e9 de qgrad y qgrad_ext frente a la prediccion registrada en el README (2.30-2.31 y 2.50-2.54 GPa) y la sensibilidad m a la tasa de deformacion.
Barras = desvio entre x, y, z del MISMO vidrio. Uso (venv_structural): python trends_qgrad.py"""
import os, csv
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SAMPLES = ["2e10", "1e11", "1e12", "qgrad", "qgrad_ext"]
COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad": "#7b4fc4", "qgrad_ext": "#b5338a"}
MK = {"2e10": "o", "1e11": "s", "1e12": "^", "qgrad": "D", "qgrad_ext": "v"}; LS = {"qgrad": "--", "qgrad_ext": "--"}
PAN = [("E_GPa", "E (GPa)"), ("sigma_max_GPa", "σmax (GPa)"), ("sigma_flujo_GPa", "σ flujo, ε 0.15–0.20 (GPa)"), ("delta_sigma_GPa", "Δσ (GPa)"), ("psi_w02_04", "ψ(0.02–0.04)"), ("eps_sigma_max", "ε(σmax)")]
g = defaultdict(list)
for r in csv.DictReader(open(f"{M}/resumen.csv")):
    if float(r["T"]) == 300.0 and r["sample"] in SAMPLES: g[(r["sample"], r["erate"])].append(r)
ms = lambda rs, c: (lambda v: (v.mean(), v.std(ddof=1) if len(v) > 1 else float("nan")))(np.array([float(r[c]) for r in rs]))
with open(f"{M}/tendencias_con_qgrad.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample", "T", "erate", "n_dirs"] + [x for c, _ in PAN for x in (c + "_media", c + "_sd")])
    for k in sorted(g, key=lambda k: (SAMPLES.index(k[0]), float(k[1]))):
        row = [k[0], 300, k[1], len(g[k])]
        for c, _ in PAN: m, s = ms(g[k], c); row += [f"{m:.5g}", f"{s:.3g}"]
        w.writerow(row)
fg, axs = plt.subplots(2, 3, figsize=(14, 7.2), constrained_layout=True)
for ax, (c, nm) in zip(axs.flat, PAN):
    for s in SAMPLES:
        pts = sorted((float(k[1]), *ms(g[k], c)) for k in g if k[0] == s)
        if pts: x, m, sd = zip(*pts); ax.errorbar(x, m, yerr=np.nan_to_num(sd), color=COL[s], marker=MK[s], ms=5.5, lw=1.5, ls=LS.get(s, "-"), capsize=3, label=s)
    ax.set(xlabel="tasa de deformación (s⁻¹)", ylabel=nm, xscale="log"); ax.grid(color="#e6e5e1", lw=.7); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axs.flat[0].legend(title="enfriamiento (K/s)", frameon=False, fontsize=9)
fg.suptitle("Dependencia con la tasa de deformación (300 K), con qgrad y qgrad_ext\nBarras: desvío entre x, y, z del mismo vidrio (no incluye variación entre vidrios)", fontsize=10.5)
fg.savefig(f"{M}/figs/tendencias_erate_con_qgrad.png", dpi=140); plt.close(fg)
print("== sigma_max (GPa, media ± sd sobre x,y,z), 300 K")
for er in ("1e8", "1e9", "1e10"):
    print(f"  erate {er}: " + "  ".join(f"{s}: {ms(g[(s, er)], 'sigma_max_GPa')[0]:.3f}±{ms(g[(s, er)], 'sigma_max_GPa')[1]:.3f}" for s in SAMPLES if (s, er) in g))
print("== prediccion README a 1e9: qgrad 2.30-2.31, qgrad_ext 2.50-2.54;  observado:",
      ", ".join(f"{s} {ms(g[(s, '1e9')], 'sigma_max_GPa')[0]:.3f}" for s in ("qgrad", "qgrad_ext") if (s, "1e9") in g))
print("== m = d ln sigma_max / d ln erate")
for s in SAMPLES:
    pts = sorted((float(k[1]), ms(g[k], "sigma_max_GPa")[0]) for k in g if k[0] == s)
    if len(pts) >= 2: print(f"  {s}: m = {np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)[0]:.4f}")
for c, nm in (("E_GPa", "E"), ("sigma_flujo_GPa", "sigma_flujo"), ("delta_sigma_GPa", "delta_sigma"), ("psi_w02_04", "psi(0.02-0.04)"), ("eps_sigma_max", "eps(sigma_max)")):
    print(f"== {nm} a 1e9: " + "  ".join(f"{s} {ms(g[(s, '1e9')], c)[0]:.4g}" for s in SAMPLES if (s, "1e9") in g))
