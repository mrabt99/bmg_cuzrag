#!/usr/bin/env python3
"""Media y desvio (ddof = 1) sobre las tres direcciones x, y, z de sigma_max, delta_sigma, E, psi(0.02-0.04) y psi(0.04-0.06), por (muestra, T, erate), a partir de resumen.csv,
y comparacion entre muestras al mismo (T, erate): diferencia de medias, error estandar combinado (n = 3 por grupo, MUY POCOS datos) y el cociente |dif|/SE.
Uso: python summary_dirs.py [--T 300] [--erate 1e9] [--a 2e10 --b 1e12]"""
import sys, csv, math, itertools
from collections import defaultdict
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
Tq, Eq = arg("--T"), arg("--erate"); A, B = arg("--a", "2e10"), arg("--b", "1e12")
Q = [("sigma_max_GPa", "sigma_max [GPa]"), ("delta_sigma_GPa", "delta_sigma [GPa]"), ("E_GPa", "E [GPa]"), ("psi_w02_04", "psi(0.02-0.04)"), ("psi_w04_06", "psi(0.04-0.06)")]
g = defaultdict(list)
import os
for r in csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "resumen.csv"))):
    if (Tq is None or r["T"] == Tq) and (Eq is None or r["erate"] == Eq): g[(r["sample"], r["T"], r["erate"])].append(r)
ms = lambda v: (sum(v) / len(v), math.sqrt(sum((x - sum(v) / len(v)) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else float("nan"))
for key in sorted(g):
    n = len(g[key]); print(f"\n== {key[0]}  T={key[1]} K  erate={key[2]}  (n = {n} direcciones: {','.join(r['dir'] for r in g[key])})")
    for c, nm in Q:
        v = [float(r[c]) for r in g[key]]; m, s = ms(v); print(f"   {nm:18s} media = {m:.4f}   desvio = {s:.4f}   valores = {', '.join(f'{x:.4f}' for x in v)}")
for (a, b) in itertools.combinations(sorted({k[0] for k in g}), 2):
    for (T, e) in sorted({(k[1], k[2]) for k in g}):
        if (a, T, e) in g and (b, T, e) in g and {a, b} == {A, B}:
            print(f"\n== {a} vs {b}  T={T} K  erate={e}:  diferencia de medias (b - a), error estandar combinado, cociente")
            for c, nm in Q:
                va, vb = [float(r[c]) for r in g[(a, T, e)]], [float(r[c]) for r in g[(b, T, e)]]; (ma, sa), (mb, sb) = ms(va), ms(vb)
                se = math.sqrt(sa ** 2 / len(va) + sb ** 2 / len(vb)) if len(va) > 1 and len(vb) > 1 else float("nan")
                print(f"   {nm:18s} {mb - ma:+.4f}  SE = {se:.4f}  |dif|/SE = {abs(mb - ma) / se if se == se and se > 0 else float('nan'):.1f}")
