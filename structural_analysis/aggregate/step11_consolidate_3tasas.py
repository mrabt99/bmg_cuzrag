#!/usr/bin/env python3
"""Tarea R (paso 1) -- consolida las 3 tasas (2e10, 1e11, 1e12) en formato tidy
rate,T,especie,tipo_cluster,fraccion,N  ->  tabla_compare_rates_3tasas.csv

Tipos (solo los validados con estadistica):
  Cu  ico_boo        <- ../boo/<tasa>/NN.boo.dump, classify/boo_types.mask_ico (q6>=0.6, W6hat<=-0.16, nnn=12)
  Cu  ico_voronoi    <- ../voronoi/abundance_voronoi.csv (OVITO default) ico
  Zr  Z14/Z15/Z16_voronoi (indice exacto <0,0,12,k>, default; Z16 = limite superior, Lu et al. PRB 102, 184201)
  Ag  ico_voronoi
BOO y Voronoi NO se combinan. Sin CN== en ningun paso. Zr-ico y Ag-FK no se incluyen.
Antes de consolidar verifica forma/conteos de las 3 tasas en ambas fuentes.
"""
import sys, glob, re
import numpy as np, pandas as pd
sys.path.insert(0, "../classify")
from boo_types import load_boo_dump, mask_ico, SPECIES

RATES = ["2e10", "1e11", "1e12"]
NT = 17
NSP = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
T_OF = lambda nn: 1900 - (nn - 1) * 100
fail = []

def check(cond, msg):
    print(("  OK   " if cond else "  FALLA ") + msg)
    if not cond: fail.append(msg)

print("== manifest (17 restarts, 1900K ... glass300K) ==")
for r in RATES:
    m = [l.strip() for l in open(f"../{r}/manifest.txt") if l.strip()]
    segs = [int(re.search(r"seg(\d+)$", x).group(1)) for x in m if re.search(r"seg\d+$", x)]
    if r == "1e12":   # templado desde 2000 K (in.quench_1e12): seg1 termina a 1900 K, sin restart sweep
        ts = [int(open(f"../{r}/snaps/{nn:02d}.dump").readlines()[1]) for nn in range(1, 17)]
        check(len(m) == NT and segs == list(range(1, 17)) and "glass300K" in m[-1] and ts == [100000 * k for k in range(1, 17)],
              f"{r}: 17 entradas = seg1..seg16 + glass300K; timesteps 100000*k (1 K/ps desde 2000 K -> config 1 = 1900 K)")
    else:             # restart sweep 1900K + seg1..seg15 + glass300K
        check(len(m) == NT and m[0].endswith("1900K") and segs == list(range(1, 16)) and "glass300K" in m[-1],
              f"{r}: 17 entradas = sweep 1900K + seg1..seg15 + glass300K")

print("\n== BOO dumps: forma y conteos ==")
boo_rows, box = [], {}
for r in RATES:
    files = sorted(glob.glob(f"../boo/{r}/*.boo.dump"))
    check(len(files) == NT, f"{r}: {len(files)} dumps BOO")
    for nn in range(1, NT + 1):
        path = f"../boo/{r}/{nn:02d}.boo.dump"
        with open(path) as f:
            hdr = [next(f) for _ in range(9)]
        n_hdr = int(hdr[3]); lo, hi = map(float, hdr[5].split()[:2]); box[(r, nn)] = hi - lo
        d = load_boo_dump(path)
        ids, types = d[:, 0].astype(int), d[:, 1].astype(int)
        ok = (n_hdr == 65536 and len(d) == 65536 and np.array_equal(np.sort(ids), np.arange(1, 65537))
              and {SPECIES[t]: int((types == t).sum()) for t in SPECIES} == NSP and np.isfinite(d[:, 2:]).all())
        if not ok: check(False, f"{r} {nn:02d}: N/ids/tipos/NaN")
        ico = mask_ico(d)
        n = int((ico & (types == 2)).sum())
        boo_rows.append((r, T_OF(nn), "Cu", "ico_boo", n / NSP["Cu"], n))
print("  (solo se imprimen las fallas por dump; 51 dumps revisados)")

print("\n== caja: lado decreciente con el indice en cada tasa (T nominal coherente) ==")
for r in RATES:
    L = [box[(r, nn)] for nn in range(1, NT + 1)]
    check(all(L[i] > L[i + 1] for i in range(NT - 1)), f"{r}: L {L[0]:.2f} -> {L[-1]:.2f} A monotono decreciente")

print("\n== Voronoi (abundance_voronoi.csv) ==")
v = pd.read_csv("../voronoi/abundance_voronoi.csv", dtype={"rate": str})
for r in RATES:
    s = v[v.rate == r]
    check(len(s) == NT * 3 * 4 and s["T"].nunique() == NT and not s.fraccion.isna().any(),
          f"{r}: {len(s)} filas (esperado 204), {s['T'].nunique()} T, sin NaN")
v["N"] = (v.fraccion * v.especie.map(NSP)).round().astype(int)
check(((v.fraccion * v.especie.map(NSP)) - v.N).abs().max() <= 0.5e-6 * 31130 + 1e-9, "conteos enteros recuperables")

print("\n== cruce BOO recalculado vs tabla_clusters_51.csv (N_ico_BOO, Cu) ==")
t = pd.read_csv("tabla_clusters_51.csv", dtype={"tasa": str, "config": str})
t = t[t.especie == "Cu"].set_index(["tasa", "T_nominal"]).N_ico_BOO
b = pd.DataFrame(boo_rows, columns=["rate", "T", "especie", "tipo_cluster", "fraccion", "N"])
mism = [(r.rate, r.T, r.N, t[(r.rate, r.T)]) for r in b.itertuples() if t[(r.rate, r.T)] != r.N]
check(not mism, f"51/51 conteos Cu-ico-BOO coinciden con la tabla previa" + (f" (difieren: {mism[:5]})" if mism else ""))

if fail:
    print("\nHAY FALLAS -> no se escribe el CSV"); sys.exit(1)

sel = [("Cu", "ico"), ("Zr", "Z14"), ("Zr", "Z15"), ("Zr", "Z16"), ("Ag", "ico")]
vv = pd.concat([v[(v.especie == sp) & (v.tipo_cluster == f"{tc}_voronoi")] for sp, tc in sel])
vv = vv[["rate", "T", "especie", "tipo_cluster", "fraccion", "N"]].copy()
vv["tipo_cluster"] = vv.tipo_cluster.str.replace("^(Z14|Z15|Z16)_voronoi$", r"\1_voronoi", regex=True)
out = pd.concat([b, vv])
out["rate"] = pd.Categorical(out.rate, RATES)
out = out.sort_values(["especie", "tipo_cluster", "rate", "T"], ascending=[True, True, True, False])
out.to_csv("tabla_compare_rates_3tasas.csv", index=False, float_format="%.6f")
print(f"\n-> tabla_compare_rates_3tasas.csv: {len(out)} filas; combos:",
      sorted(set(zip(out.especie, out.tipo_cluster))))
