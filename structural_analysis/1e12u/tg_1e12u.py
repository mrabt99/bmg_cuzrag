#!/usr/bin/env python3
"""Tg del 1e12 UNIFICADO (rate_1e12/, mismo protocolo que 2e10 y 1e11) con compute_tg.py tal cual
(mismas ventanas: líquido 900-1200 K, vítreo 300-600 K, bines de 10 K) y comparación con la referencia
vieja (1e12/thermo_1e12.csv, rate_1e12_reference). No toca Tg_resultados.csv: escribe 1e12u/Tg_1e12u.csv.
Uso (desde structural_analysis/): PYTHONDONTWRITEBYTECODE=1 /scratch/negishi/aaabt/venv_structural/bin/python 1e12u/tg_1e12u.py
"""
import csv, sys
import numpy as np
sys.path.insert(0, ".")
from compute_tg import load_thermo, bin_by_temperature, two_branch_tg, TG_WINDOWS

N = 65536
out = []
for nombre, f in (("1e12 unificada", "1e12u/thermo_1e12u.csv"), ("1e12 referencia vieja", "1e12/thermo_1e12.csv"),
                  ("1e11", "1e11/thermo_1e11.csv"), ("2e10", "2e10/thermo_2e10.csv")):
    T, Ep, V = load_thermo(f)
    Tb, Vb = bin_by_temperature(T, V, 10.0); _, Epb = bin_by_temperature(T, Ep, 10.0)
    tv = two_branch_tg(Tb, Vb, TG_WINDOWS["liquid"], TG_WINDOWS["glass"])[0]
    te = two_branch_tg(Tb, Epb, TG_WINDOWS["liquid"], TG_WINDOWS["glass"])[0]
    d = list(csv.DictReader(open(f)))
    fin = d[-50:]   # últimas 50 filas = final de la relajación a 300 K
    ep = np.mean([float(r["Ep"]) for r in fin]); rho = np.mean([float(r["rho"]) for r in fin])
    lx = np.mean([float(r["Lx"]) for r in fin]); tt = np.mean([float(r["T"]) for r in fin])
    out.append([nombre, tv, te, (tv + te) / 2, tt, ep, rho, lx])
    print(f"{nombre:24s} Tg(V)={tv:6.1f}  Tg(Ep)={te:6.1f}  media={(tv+te)/2:6.1f} K | final: T={tt:6.1f} K "
          f"Ep={ep:.5f} eV  rho={rho:.4f} g/cc  Lx={lx:.3f} A")
with open("1e12u/Tg_1e12u.csv", "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["muestra", "Tg_V_K", "Tg_Ep_K", "Tg_media_K", "T_final_K", "Ep_final_eV", "rho_final", "Lx_final"])
    for r in out:
        w.writerow([r[0]] + [f"{x:.1f}" for x in r[1:5]] + [f"{r[5]:.5f}", f"{r[6]:.4f}", f"{r[7]:.3f}"])
