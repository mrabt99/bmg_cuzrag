#!/usr/bin/env python3
"""Analisis de la curva sigma-epsilon de un caso de la matriz (version de pilot_tension/analyze_curve.py generalizada al eje de traccion).
Uso: python analyze_curve.py <caso>/<tag>.curve.dat --dir x|y|z [--out <carpeta>] [--label texto]
Columnas del .curve.dat: step strain sxx syy szz T PE lx ly lz vol (GPa; traccion > 0). La tension del eje de traccion es sxx, syy o szz segun --dir.
Reporta E (ajuste lineal strain <= 0.01, ordenada libre), sigma_max y su strain, sigma_flujo (media 0.15-0.20), delta_sigma = sigma_max - sigma_flujo.
Salidas en <out> (por defecto <carpeta del curve>/analysis): curva_resumen.txt, curva_resumen.json, curva_sigma_eps.png.
Los errores son SOLO ruido termico de una realizacion."""
import sys, os, json
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
path = sys.argv[1]; DIR = arg("--dir", "x"); out = arg("--out", os.path.join(os.path.dirname(os.path.abspath(path)), "analysis")); label = arg("--label", os.path.basename(os.path.dirname(os.path.abspath(path))))
os.makedirs(out, exist_ok=True); ci = {"x": 2, "y": 3, "z": 4}[DIR]
d = np.loadtxt(path, comments="#"); e, s = d[:, 1], d[:, ci]; T = d[:, 5]
lat = [c for c in (2, 3, 4) if c != ci]
EMAX_E, F0, F1, NSM = 0.01, 0.15, 0.20, 5
m = e <= EMAX_E; (a, b), cov = np.polyfit(e[m], s[m], 1, cov=True); a_err = np.sqrt(cov[0, 0])
i = int(np.argmax(s)); sm = np.convolve(s, np.ones(NSM) / NSM, mode="valid"); es = np.convolve(e, np.ones(NSM) / NSM, mode="valid"); j = int(np.argmax(sm))
f = (e >= F0) & (e <= F1); sf, sf_sd = s[f].mean(), s[f].std(ddof=1)
res = dict(label=label, dir=DIR, n=int(len(e)), E=float(a), E_err=float(a_err), sigma_max=float(s[i]), eps_sigma_max=float(e[i]), sigma_max_smooth=float(sm[j]), eps_sigma_max_smooth=float(es[j]),
           sigma_flow=float(sf), sigma_flow_sd=float(sf_sd), delta_sigma=float(s[i] - sf), lateral_max=float(max(np.abs(d[:, lat[0]]).max(), np.abs(d[:, lat[1]]).max())), T_mean=float(T.mean()), T_sd=float(T.std()))
lines = [f"{label}: eje {DIR} ; archivo {os.path.basename(path)} ({len(e)} puntos, strain {e[0]:.4f} - {e[-1]:.4f})",
         f"E (ajuste lineal, strain <= {EMAX_E}, {int(m.sum())} puntos) = {a:.1f} +- {a_err:.1f} GPa",
         f"sigma_max = {s[i]:.3f} GPa en strain = {e[i]:.4f}   (suavizado 5 puntos: {sm[j]:.3f} GPa en {es[j]:.4f})",
         f"sigma_flujo (media {F0}-{F1}, {int(f.sum())} puntos) = {sf:.3f} +- {sf_sd:.3f} GPa ; delta_sigma = {s[i]-sf:.3f} GPa",
         f"verificacion: |tensiones laterales| max = {res['lateral_max']:.4f} GPa ; T = {T.mean():.2f} +- {T.std():.2f} K"]
print("\n".join(lines)); open(f"{out}/curva_resumen.txt", "w").write("\n".join(lines) + "\n"); json.dump(res, open(f"{out}/curva_resumen.json", "w"), indent=1)
BL, GR, INK = "#2a78d6", "#8a8985", "#0b0b0b"
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw=dict(width_ratios=[2.2, 1])); A = ax[0]
A.plot(e, s, color=BL, lw=1.6, label=f"σ{DIR}{DIR} (promedio por Δε = 0.001)"); ee = np.array([0, 0.012]); A.plot(ee, a * ee + b, color=GR, ls="--", lw=1.2, label=f"ajuste E = {a:.0f} GPa (ε ≤ {EMAX_E})")
A.hlines(sf, F0, F1, color=INK, lw=2.0, label=f"σ flujo = {sf:.2f} GPa (ε {F0}–{F1})"); A.plot([e[i]], [s[i]], "o", color=INK, ms=6, zorder=5)
A.annotate(f"σmax = {s[i]:.2f} GPa\nε = {e[i]:.3f}", (e[i], s[i]), xytext=(e[i] + 0.012, s[i] + 0.02), fontsize=9, color=INK)
A.set(xlabel="deformación ingeniería ε", ylabel=f"σ{DIR}{DIR} (GPa)", xlim=(0, 0.205), ylim=(0, max(s.max(), 0.1) * 1.15), title=label); A.grid(color="#e6e5e1", lw=.7); A.set_axisbelow(True); A.legend(frameon=False, fontsize=8.5, loc="lower right")
B = ax[1]; B.plot(e, d[:, lat[0]], color="#eb6834", lw=1.1, label="σ lateral 1"); B.plot(e, d[:, lat[1]], color="#1baf7a", lw=1.1, label="σ lateral 2"); B.axhline(0, color=GR, lw=.8)
B.set(xlabel="ε", ylabel="σ (GPa)", title="Verificación: tensiones laterales (≈ 0)", xlim=(0, 0.205)); B.grid(color="#e6e5e1", lw=.7); B.legend(frameon=False, fontsize=8.5)
for x in ax:
    for sp in ("top", "right"): x.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(f"{out}/curva_sigma_eps.png", dpi=150)
