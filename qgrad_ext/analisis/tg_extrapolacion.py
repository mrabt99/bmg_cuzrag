#!/usr/bin/env python3
"""Paso 1: Tg vs log10(tasa) de las 3 tasas (Etapa 1) y extrapolacion a 1.5e9 K/s; T_switch = Tg - 100 K (o -50 K) redondeado a multiplo de 50 K.
Lee (solo lectura) structural_analysis/qgrad/thermo/Tg_qgrad_y_sensibilidad.csv, que reproduce exacto los Tg de compute_tg.py (646.1/637.2, 669.3/663.0, 703.2/697.9).
SUPUESTOS: solo 3 puntos; la extrapolacion (1.125 decadas por debajo de la tasa mas lenta) asume linealidad en log10(tasa). El error estandar del ajuste (1 g.l.) NO es la
incertidumbre real; la real es la sensibilidad a las ventanas de ajuste (se reporta el rango). Uso: PYTHONDONTWRITEBYTECODE=1 python3 tg_extrapolacion.py"""
import csv, numpy as np
rows = [r for r in csv.DictReader(open("/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis/qgrad/thermo/Tg_qgrad_y_sensibilidad.csv")) if r["curva"] in ("2e10", "1e11", "1e12")]
x = np.array([np.log10(2e10), 11.0, 12.0]); x0 = np.log10(1.5e9); wins = list(dict.fromkeys(r["ventanas"] for r in rows)); out = []
for w in wins:
    t = {k: next(r for r in rows if r["curva"] == k and r["ventanas"] == w) for k in ("2e10", "1e11", "1e12")}
    for nm, col in (("Tg(V)", "Tg_V_K"), ("Tg(Ep)", "Tg_Ep_K"), ("media", None)):
        y = np.array([float(t[k]["Tg_V_K"]) if col == "Tg_V_K" else float(t[k]["Tg_Ep_K"]) if col else (float(t[k]["Tg_V_K"]) + float(t[k]["Tg_Ep_K"])) / 2 for k in t])
        b, a = np.polyfit(x, y, 1); res = y - (a + b * x); s = np.sqrt((res ** 2).sum()); Sxx = ((x - x.mean()) ** 2).sum()
        se = s * np.sqrt(1 + 1 / 3 + (x0 - x.mean()) ** 2 / Sxx); b2 = (y[1] - y[0]) / (x[1] - x[0]); out.append([w, nm, *y, b, a + b * x0, se, y[0] + b2 * (x0 - x[0])])
with open("tg_extrapolacion.csv", "w", newline="") as f:
    w_ = csv.writer(f); w_.writerow(["ventanas", "tipo", "Tg_2e10", "Tg_1e11", "Tg_1e12", "pendiente_K_por_decada", "Tg_1.5e9_ajuste3", "SE_pred_1gl", "Tg_1.5e9_solo_2_lentas"]); w_.writerows([[o[0], o[1]] + [f"{v:.2f}" for v in o[2:]] for o in out])
std = [o for o in out if o[0].startswith("ventanas estandar") and o[1] == "media"][0]; med = [o[6] for o in out if o[1] == "media"]
print(f"log10(1.5e9) = {x0:.4f}\nCRITERIO Etapa 1 (ventanas estandar, media de Tg(V) y Tg(Ep)): Tg(1.5e9) = {std[6]:.1f} K (pendiente {std[5]:.1f} K/decada; SE 1 g.l. +-{std[7]:.1f} K; solo 2e10+1e11: {std[8]:.1f} K)")
print(f"Rango por eleccion de ventanas (media): {min(med):.0f}-{max(med):.0f} K")
for d in (100, 50): print(f"T_switch = Tg - {d} K = {std[6]-d:.1f} K -> multiplo de 50: {int(round((std[6]-d)/50)*50)} K")
