#!/usr/bin/env python3
"""Razon qgrad_ext / referencia en el vidrio y tabla por T. Criterios identicos a qgrad/aggregate/regimen_qgrad.py:
filtro N>=20 por T en AMBAS curvas, razon = N_pool(ext)/N_pool(ref), error = razon*sqrt(1/N1+1/N2) (piso Poisson), 300 K = 17r (relajado).
Ventana vitrea FIJA = {500, 400, 300} K (la MISMA que dio la razon de qgrad: comparacion directa; el Tg de la rama no se midio y la
ventana Tg-100 de qgrad era T<531.9 K). Referencias: 2e10, 1e11, 1e12, qgrad.
Salidas: tabla_regimen_ext.csv (pool 500/400/300) y tabla_porT_ext.csv (N y fraccion de 2e10, qgrad, qgrad_ext a 700..300 K + razones).
Uso: python regimen_ext.py"""
import csv, math
NMIN, WIN = 20, (500, 400, 300)
N = {}
for r in csv.DictReader(open("tabla_compare_rates_5curvas.csv")):
    T = int(r["T"]); rate = r["rate"]
    if T == 300 and rate in ("qgrad", "qgrad_ext"):
        if r["config"] != "17r": continue
    N[(rate, T, r["especie"], r["tipo_cluster"])] = int(r["N"])
combos = sorted({(k[2], k[3]) for k in N})
out = []
for ref in ("2e10", "1e11", "1e12", "qgrad"):
    for sp, tc in combos:
        per = [(T, N.get(("qgrad_ext", T, sp, tc)), N.get((ref, T, sp, tc))) for T in WIN]
        per = [(T, a, b) for T, a, b in per if a is not None and b is not None]
        kept = [(T, a, b) for T, a, b in per if a >= NMIN and b >= NMIN]
        n1, n2 = sum(a for _, a, _ in kept), sum(b for _, _, b in kept)
        if n1 and n2:
            q = n1 / n2; e = q * math.sqrt(1 / n1 + 1 / n2); v = (n1, n2, f"{q:.4f}", f"{e:.4f}", f"{(q - 1) / e:.2f}")
        else: v = ("", "", "", "", "")
        out.append((ref, sp, tc, " ".join(str(t) for t, _, _ in kept), *v, "" if kept else "NO CALCULABLE (N<20)"))
with open("tabla_regimen_ext.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["referencia", "especie", "tipo_cluster", "T_usados", "N_pool_ext", "N_pool_ref", "razon_ext/ref", "error", "(razon-1)/error", "nota"]); w.writerows(out)
KEY = [("Cu", "ico_boo"), ("Cu", "ico_voronoi"), ("Zr", "Z14_voronoi"), ("Zr", "Z15_voronoi"), ("Zr", "Z16_voronoi"), ("Ag", "ico_voronoi")]
print("=== razon vitrea (500/400/300 K, 300=17r) qgrad_ext / ref ===")
for ref in ("2e10", "qgrad", "1e11", "1e12"):
    print(f"-- ref = {ref}")
    for o in out:
        if o[0] == ref and (o[1], o[2]) in KEY and o[6]: print(f"   {o[1]:2s} {o[2]:12s} T[{o[3]}] pool {o[4]}/{o[5]}  razon {o[6]} +- {o[7]}  ({o[8]} sigma)")
with open("tabla_porT_ext.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["especie", "tipo_cluster", "T", "N_2e10", "N_qgrad", "N_ext", "ext/2e10", "ext/qgrad", "qgrad/2e10"])
    print("\n=== por T: N (2e10 / qgrad / ext) y razones ===")
    for sp, tc in KEY:
        print(f"{sp} {tc}")
        for T in (800, 700, 600, 500, 400, 300):
            a, b, c = (N.get((r, T, sp, tc)) for r in ("2e10", "qgrad", "qgrad_ext"))
            if None in (a, b, c): continue
            f1 = lambda x, y: f"{x / y:.3f}" if y else ""
            w.writerow([sp, tc, T, a, b, c, f1(c, a), f1(c, b), f1(b, a)])
            print(f"   T={T}: {a:5d} {b:5d} {c:5d}   ext/2e10 {f1(c, a)}  ext/qgrad {f1(c, b)}  qgrad/2e10 {f1(b, a)}")
