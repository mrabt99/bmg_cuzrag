#!/usr/bin/env python3
"""Analisis de la curva sigma-epsilon del ensayo de traccion (lee <tag>.curve.dat de LAMMPS).
Columnas del .curve.dat: step strain sxx syy szz T PE(eV/at) lx ly lz vol  (tensiones en GPa, traccion > 0).
Reporta:
  E        pendiente del ajuste lineal sxx(strain) para strain <= EMAX_E (0.01); con ordenada libre (principal) y por el origen
  sigma_max y su strain (maximo de la curva promediada; tambien del suavizado con media movil de 5 puntos)
  sigma_flujo  media de sxx para 0.15 <= strain <= 0.20 (y su desvio)
  dsigma   = sigma_max - sigma_flujo
Salidas (en analysis/): curva_sigma_eps.png, curva_resumen.txt.  Uso: python analyze_curve.py [archivo.curve.dat]
Los errores (error estandar del ajuste, desvio dentro de la ventana) son solo ruido termico de UNA realizacion: no incluyen variacion entre muestras."""
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
path = sys.argv[1] if len(sys.argv) > 1 else "tens_2e10_300K_1e9_x.curve.dat"
os.makedirs("analysis", exist_ok=True)
d = np.loadtxt(path, comments="#"); e, sxx, syy, szz, T = d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
EMAX_E, F0, F1, NSM = 0.01, 0.15, 0.20, 5

m = e <= EMAX_E
(a, b), cov = np.polyfit(e[m], sxx[m], 1, cov=True); a_err = np.sqrt(cov[0, 0])
E0 = np.sum(e[m] * sxx[m]) / np.sum(e[m] ** 2)                     # ajuste por el origen
i = int(np.argmax(sxx)); sm = np.convolve(sxx, np.ones(NSM) / NSM, mode="valid"); es = np.convolve(e, np.ones(NSM) / NSM, mode="valid"); j = int(np.argmax(sm))
f = (e >= F0) & (e <= F1); sf, sf_sd = sxx[f].mean(), sxx[f].std(ddof=1)
lines = [f"archivo: {path}   ({len(e)} puntos, strain {e[0]:.4f} - {e[-1]:.4f})",
         f"E (ajuste lineal, strain <= {EMAX_E}, {int(m.sum())} puntos, ordenada libre) = {a:.1f} +- {a_err:.1f} GPa   (ordenada {b:.3f} GPa)",
         f"E (ajuste por el origen, mismos puntos)                                  = {E0:.1f} GPa",
         f"sigma_max (curva promediada, maximo)                  = {sxx[i]:.3f} GPa  en strain = {e[i]:.4f}",
         f"sigma_max (suavizado, media movil de {NSM} puntos)        = {sm[j]:.3f} GPa  en strain = {es[j]:.4f}",
         f"sigma_flujo (media, {F0} <= strain <= {F1}, {int(f.sum())} puntos)  = {sf:.3f} +- {sf_sd:.3f} GPa (desvio dentro de la ventana)",
         f"delta_sigma = sigma_max - sigma_flujo                 = {sxx[i] - sf:.3f} GPa   (con sigma_max suavizado: {sm[j] - sf:.3f})",
         f"verificacion: |syy|, |szz| maximos = {np.abs(syy).max():.4f}, {np.abs(szz).max():.4f} GPa ; media = {syy.mean():+.4f}, {szz.mean():+.4f} GPa ; T = {T.mean():.2f} +- {T.std():.2f} K (min {T.min():.1f}, max {T.max():.1f})",
         f"caja: lx {d[0,7]:.2f} -> {d[-1,7]:.2f} ; ly {d[0,8]:.2f} -> {d[-1,8]:.2f} ; lz {d[0,9]:.2f} -> {d[-1,9]:.2f} A ; vol {d[0,10]:.0f} -> {d[-1,10]:.0f} A^3 ({100*(d[-1,10]/d[0,10]-1):+.2f} %)"]
print("\n".join(lines)); open("analysis/curva_resumen.txt", "w").write("\n".join(lines) + "\n")

BL, GR, INK = "#2a78d6", "#8a8985", "#0b0b0b"
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw=dict(width_ratios=[2.2, 1]))
A = ax[0]
A.plot(e, sxx, color=BL, lw=1.6, label="σxx (promedio por Δε = 0.001)")
ee = np.array([0, 0.012]); A.plot(ee, a * ee + b, color=GR, ls="--", lw=1.2, label=f"ajuste E = {a:.0f} GPa (ε ≤ {EMAX_E})")
A.hlines(sf, F0, F1, color=INK, lw=2.0, label=f"σ flujo = {sf:.2f} GPa (ε {F0}–{F1})")
A.plot([e[i]], [sxx[i]], "o", color=INK, ms=6, zorder=5); A.annotate(f"σmax = {sxx[i]:.2f} GPa\nε = {e[i]:.3f}", (e[i], sxx[i]), xytext=(e[i] + 0.012, sxx[i] + 0.02), fontsize=9, color=INK)
A.set(xlabel="deformación ingeniería ε", ylabel="σxx (GPa)", xlim=(0, 0.205), ylim=(0, sxx.max() * 1.15), title="Tracción uniaxial en x — 2e10, 300 K, 1e9 s⁻¹")
A.grid(color="#e6e5e1", lw=.7); A.set_axisbelow(True); A.legend(frameon=False, fontsize=8.5, loc="lower right")
B = ax[1]; B.plot(e, syy, color="#eb6834", lw=1.1, label="σyy"); B.plot(e, szz, color="#1baf7a", lw=1.1, label="σzz"); B.axhline(0, color=GR, lw=.8)
B.set(xlabel="ε", ylabel="σ (GPa)", title="Verificación: tensiones laterales (≈ 0)", xlim=(0, 0.205)); B.grid(color="#e6e5e1", lw=.7); B.legend(frameon=False, fontsize=8.5)
for x in ax:
    for s in ("top", "right"): x.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("analysis/curva_sigma_eps.png", dpi=150); print("-> analysis/curva_sigma_eps.png, analysis/curva_resumen.txt")
