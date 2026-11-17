#!/usr/bin/env python3
"""Tg del PERFIL qgrad_ext (2e10->1.5e9 K/s en 1500-800 K, 1.5e9 en 800-500 K, 1e11 en 500-300 K) con la MISMA logica y
ventanas que compute_tg.py / qgrad/thermo/tg_qgrad.py (doble ajuste lineal, bines de 10 K; estandar: liq 900-1200, vit 300-600)
y las 3 ventanas alternativas de sensibilidad. Luego compara con Tg(1.5e9) extrapolada (tg_extrapolacion.csv) e invierte la recta
Tg(log10 gamma) de las 3 tasas constantes para dar una "tasa efectiva por Tg" (misma ventana en ambos lados).

Entrada: thermo_ext.csv (extract_thermo.py --style avetime sobre qgrad seg01-08 + ext seg09,10a,10,11,12; las filas de seg01-08
son identicas a structural_analysis/qgrad/thermo/thermo_qgrad.csv, verificado). Solo LEE fuera de esta carpeta.
AVISO: la ventana vitrea estandar 300-600 K mezcla 600-500 K (1.5e9) con 500-300 K (1e11); la pendiente vitrea casi no depende
de la tasa (0.143-0.149 meV/at/K en las 5 curvas), pero es un supuesto. Tg del perfil != Tg de una tasa constante.
Uso: PYTHONDONTWRITEBYTECODE=1 /scratch/negishi/aaabt/venv_structural/bin/python tg_ext.py
"""
import sys, csv
import numpy as np
SA = "/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis"
sys.path.insert(0, SA)
from compute_tg import load_thermo, bin_by_temperature, two_branch_tg, TG_WINDOWS

ALT = {"ventanas estandar (liq 900-1200 / vit 300-600)": TG_WINDOWS,
       "liq 900-1100 / vit 300-500": {"liquid": (900., 1100.), "glass": (300., 500.)},
       "liq 800-1000 / vit 300-600": {"liquid": (800., 1000.), "glass": (300., 600.)},
       "liq 1000-1200 / vit 300-600": {"liquid": (1000., 1200.), "glass": (300., 600.)}}

T, Ep, V = load_thermo("thermo_ext.csv")
Tb, Vb = bin_by_temperature(T, V, 10.0); _, Epb = bin_by_temperature(T, Ep, 10.0)
ref = {r["ventanas"]: r for r in csv.DictReader(open(f"{SA}/qgrad/thermo/Tg_qgrad_y_sensibilidad.csv")) if r["curva"] == "qgrad"}
ext = {r["ventanas"]: r for r in csv.DictReader(open("tg_extrapolacion.csv")) if r["tipo"] == "media"}
x = np.array([np.log10(2e10), 11.0, 12.0])

rows = []
for name, w in ALT.items():
    tv = two_branch_tg(Tb, Vb, w["liquid"], w["glass"])[0]
    te = two_branch_tg(Tb, Epb, w["liquid"], w["glass"])[0]
    tm = (tv + te) / 2
    e = ext[name]
    y = np.array([float(e["Tg_2e10"]), float(e["Tg_1e11"]), float(e["Tg_1e12"])])
    b, a = np.polyfit(x, y, 1)
    g_eff = 10 ** ((tm - a) / b)
    q = ref[name]; tq = (float(q["Tg_V_K"]) + float(q["Tg_Ep_K"])) / 2
    gq = 10 ** ((tq - a) / b)
    rows.append([name, tv, te, tm, float(e["Tg_1.5e9_ajuste3"]), tm - float(e["Tg_1.5e9_ajuste3"]), g_eff, tq, gq])

with open("Tg_ext_y_sensibilidad.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["ventanas", "Tg_V_ext", "Tg_Ep_ext", "Tg_media_ext", "Tg_1.5e9_extrapolada", "ext_menos_extrap",
                 "gamma_eff_por_Tg_ext", "Tg_media_qgrad", "gamma_eff_por_Tg_qgrad"])
    for r in rows:
        wr.writerow([r[0]] + [f"{v:.1f}" for v in r[1:6]] + [f"{r[6]:.3g}", f"{r[7]:.1f}", f"{r[8]:.3g}"])
for r in rows:
    print(f"{r[0]:48s} ext Tg(V)={r[1]:6.1f} Tg(Ep)={r[2]:6.1f} media={r[3]:6.1f} | extrap(1.5e9)={r[4]:6.1f} "
          f"diff={r[5]:+5.1f} | g_eff={r[6]:.2e} | qgrad media={r[7]:6.1f} g_eff={r[8]:.2e}")
