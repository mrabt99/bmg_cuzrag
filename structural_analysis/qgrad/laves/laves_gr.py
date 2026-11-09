#!/usr/bin/env python3
"""g(r) parciales (Zr-Zr, Cu-Cu, Zr-Cu; run 0 de un unico snapshot por T): diferencia de qgrad respecto de 2e10 (control) comparada con la
diferencia 1e12-2e10 y 1e11-2e10 (escala 'normal' de cambio por tasa) y con el ruido. Un grano Laves daria picos NUEVOS y agudos
(C15 ideal a 300 K: Zr-Zr 3.07, Cu-Cu 2.51, Zr-Cu 2.94 A, luego a/sqrt2=5.0). Se reporta rms y max|dg| en r en [2,8] A.
Lee etapa4 4curvas (solo lectura). Salida laves_gr.csv"""
import csv
import numpy as np
E4 = "../etapa4/"
def load(path, cols):
    d = {}
    for r in csv.DictReader(open(path)):
        d.setdefault((r["rate"], int(float(r["T"]))), []).append([float(r["r"])] + [float(r[c]) for c in cols])
    return {k: np.array(v) for k, v in d.items()}
S = load(E4 + "rdf_sametype/rdf_sametype_all.csv", ["g_ZrZr", "g_CuCu"]); X = load(E4 + "rdf_cross/rdf_cross_all.csv", ["g_ZrCu"])
rows = []
print(f"{'par':6s} {'T':>5s}  " + "  ".join(f"{k:>18s}" for k in ("qgrad-2e10", "1e11-2e10", "1e12-2e10")) + "   (rms | max|dg|)")
for pair, D, j in (("Zr-Zr", S, 1), ("Cu-Cu", S, 2), ("Zr-Cu", X, 1)):
    for T in (700, 500, 300):
        r = D[("2e10", T)][:, 0]; m = (r >= 2) & (r <= 8); base = D[("2e10", T)][:, j]
        out = []
        for k in ("qgrad", "1e11", "1e12"):
            dg = (D[(k, T)][:, j] - base)[m]; out.append((np.sqrt((dg ** 2).mean()), np.abs(dg).max(), r[m][np.argmax(np.abs(dg))]))
        rows.append([pair, T] + [f"{v:.4f}" for o in out for v in o])
        print(f"{pair:6s} {T:5d}  " + "  ".join(f"{o[0]:.3f} | {o[1]:.3f}@{o[2]:.2f}" for o in out))
with open("laves_gr.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["par", "T"] + [f"{k}_{q}" for k in ("qgrad", "1e11", "1e12") for q in ("rms", "maxabs", "r_at_max")]); w.writerows(rows)
