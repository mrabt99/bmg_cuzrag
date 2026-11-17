#!/usr/bin/env python3
"""Tasas efectivas de qgrad y qgrad_ext recalculadas con la 1e12 UNIFICADA (rate_1e12/, etiqueta 1e12u; 17/11/2026).

Reemplaza el punto 1e12 de las dos leyes contra log10(gamma) por el de la corrida unificada:
  (1) Tg(log gamma), con las MISMAS 4 ventanas de tg_extrapolacion.py / tg_ext.py (compute_tg.py sin cambios);
      2e10, 1e11, qgrad: de structural_analysis/qgrad/thermo/Tg_qgrad_y_sensibilidad.csv; qgrad_ext: de Tg_ext_y_sensibilidad.csv;
      1e12u: calculado acá de structural_analysis/1e12u/thermo_1e12u.csv.
  (2) Ep(300 K): promedio de los últimos 5 ns de la relajación de 10 ns, para las 5 muestras. La 1e12 vieja no tenía
      log de la relajación y se usaba Ep del "run 0"; la 1e12u sí lo tiene, así que ahora el criterio es el mismo para todas.
Solo lee; escribe tasa_efectiva_1e12u.csv. Los archivos viejos (tg_extrapolacion.csv, Tg_ext_y_sensibilidad.csv) no se tocan.
Uso: PYTHONDONTWRITEBYTECODE=1 /scratch/negishi/aaabt/venv_structural/bin/python tasa_efectiva_1e12u.py
"""
import csv, sys
import numpy as np
B = "/scratch/negishi/aaabt/bmg_cuzrag"
sys.path.insert(0, f"{B}/structural_analysis")
from compute_tg import load_thermo, bin_by_temperature, two_branch_tg, TG_WINDOWS  # noqa: E402

VENT = {"ventanas estandar (liq 900-1200 / vit 300-600)": TG_WINDOWS,
        "liq 900-1100 / vit 300-500": {"liquid": (900., 1100.), "glass": (300., 500.)},
        "liq 800-1000 / vit 300-600": {"liquid": (800., 1000.), "glass": (300., 600.)},
        "liq 1000-1200 / vit 300-600": {"liquid": (1000., 1200.), "glass": (300., 600.)}}
x = np.array([np.log10(2e10), 11.0, 12.0]); x0 = np.log10(1.5e9)

# ---------------------------------------------------------------- Tg
T, Ep, V = load_thermo(f"{B}/structural_analysis/1e12u/thermo_1e12u.csv")
Tb, Vb = bin_by_temperature(T, V, 10.0); _, Epb = bin_by_temperature(T, Ep, 10.0)
tq = {}
for r in csv.DictReader(open(f"{B}/structural_analysis/qgrad/thermo/Tg_qgrad_y_sensibilidad.csv")):
    tq[(r["curva"], r["ventanas"])] = (float(r["Tg_V_K"]) + float(r["Tg_Ep_K"])) / 2
te = {r["ventanas"]: float(r["Tg_media_ext"]) for r in csv.DictReader(open("Tg_ext_y_sensibilidad.csv"))}
filas = []
print("Tg (media de V y Ep) y tasa efectiva por Tg, con la 1e12 unificada:")
for w, win in VENT.items():
    t12 = (two_branch_tg(Tb, Vb, win["liquid"], win["glass"])[0] + two_branch_tg(Tb, Epb, win["liquid"], win["glass"])[0]) / 2
    y = np.array([tq[("2e10", w)], tq[("1e11", w)], t12])
    b, a = np.polyfit(x, y, 1)
    ge = {k: 10 ** ((v - a) / b) for k, v in (("qgrad", tq[("qgrad", w)]), ("qgrad_ext", te[w]))}
    vieja = tq[("1e12", w)]
    print(f"  {w:46s} Tg(1e12u) = {t12:6.1f} (vieja {vieja:6.1f}) | pendiente {b:5.1f} K/déc | Tg(1,5e9) = {a + b * x0:6.1f} K"
          f" | γ_eff qgrad {ge['qgrad']:.2e}, ext {ge['qgrad_ext']:.2e}")
    filas.append(["Tg", w, f"{t12:.1f}", f"{vieja:.1f}", f"{b:.2f}", f"{a + b * x0:.1f}", f"{ge['qgrad']:.3g}", f"{ge['qgrad_ext']:.3g}"])

# ---------------------------------------------------------------- Ep(300 K), últimos 5 ns
RLX = {"2e10": "rate_2e10/prod_2e10/relax300K_N65536_2e10.log.dat", "1e11": "rate_1e11/prod_1e11/relax300K_N65536_1e11.log.dat",
       "1e12u": "rate_1e12/prod_1e12/relax300K_N65536_1e12.log.dat",
       "qgrad": "structural_analysis/qgrad/relax300K/relax300K_N65536_qgrad.log.dat", "qgrad_ext": "qgrad_ext/relax300K_N65536_gradext.log.dat"}
EP5 = {}
for k, f in RLX.items():
    d = np.loadtxt(f"{B}/{f}", comments="#")
    sel = d[:, 0] >= d[-1, 0] - 5_000_000
    EP5[k] = d[sel, 2].mean()
    print(f"  Ep5 {k:9s} = {EP5[k]:.6f} eV/at  (T media {d[sel, 1].mean():.1f} K, {sel.sum()} filas)")
y = np.array([EP5["2e10"], EP5["1e11"], EP5["1e12u"]])
b, a = np.polyfit(x, y, 1)
ge = {k: 10 ** ((EP5[k] - a) / b) for k in ("qgrad", "qgrad_ext")}
print(f"Ley Ep: {1e3 * b:.2f} meV/átomo por década (antes 4,32 con la 1e12 vieja run 0); "
      f"Ep predicha a 1,5e9 = {a + b * x0:.6f} eV; ΔEp ext - 2e10 = {1e3 * (EP5['qgrad_ext'] - EP5['2e10']):+.2f} meV/at")
print(f"γ_eff por Ep: qgrad {ge['qgrad']:.2e}, qgrad_ext {ge['qgrad_ext']:.2e}  (antes 1,41e10 y 2,84e9)")
print(f"faltante de ext respecto de la ley en 1,5e9: {1e3 * (EP5['qgrad_ext'] - (a + b * x0)):+.2f} meV/at")
filas.append(["Ep5", "últimos 5 ns de la relajación", f"{EP5['1e12u']:.6f}", "", f"{1e3 * b:.3f} meV/déc",
              f"{a + b * x0:.6f}", f"{ge['qgrad']:.3g}", f"{ge['qgrad_ext']:.3g}"])
with open("tasa_efectiva_1e12u.csv", "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["ley", "ventanas", "valor_1e12u", "valor_1e12_vieja", "pendiente", "pred_1.5e9", "gamma_eff_qgrad", "gamma_eff_qgrad_ext"])
    w.writerows(filas)
