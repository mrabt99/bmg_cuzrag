#!/usr/bin/env python3
"""Tarea P -- unifica la fuente Voronoi de tabla_regimen_pool.csv: ico_voronoi DEFAULT
(../voronoi/abundance_voronoi.csv, OVITO sin threshold) para Zr liq, Ag liq, Cu-Voronoi
liq+vitreo. Cu-BOO no se toca. Filtro N>=20 POR T en ambas tasas; T que falla se
descarta; si ninguno pasa -> "no calculable". Error razon*sqrt(1/N1+1/N2) sobre pool.
Ventanas compartidas: liquido = intersec(T>Tg+100), vitreo = intersec(T<Tg-100).
Uso: python step9_unify_ico_default.py [--write]
"""
import sys, math
import pandas as pd

WRITE = "--write" in sys.argv
NTOT = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
NMIN = 20
TG = {"2e10": (646.1 + 637.2) / 2, "1e12": (703.2 + 697.9) / 2}
BUF = 100.0

a = pd.read_csv("../voronoi/abundance_voronoi.csv", dtype={"rate": str})
a["tipo_cluster"] = a.tipo_cluster.str.removesuffix("_voronoi")
a["N"] = (a.fraccion * a.especie.map(NTOT)).round().astype(int)
Tall = sorted(a["T"].unique(), reverse=True)
WIN = {}
for reg, f in [("liquido", lambda r, t: t > TG[r] + BUF), ("vitreo", lambda r, t: t < TG[r] - BUF)]:
    WIN[reg] = sorted(set(t for t in Tall if f("2e10", t)) & set(t for t in Tall if f("1e12", t)), reverse=True)
print("ventanas compartidas:", {k: (v[0], v[-1], len(v)) for k, v in WIN.items()})
assert WIN["vitreo"] == [500, 400, 300]


def pooled(sp, reg):
    s = a[(a.especie == sp) & (a.tipo_cluster == "ico") & a["T"].isin(WIN[reg]) & a.rate.isin(["2e10", "1e12"])]
    p = s.pivot(index="T", columns="rate", values="N").loc[WIN[reg]]
    p["pasa"] = (p["2e10"] >= NMIN) & (p["1e12"] >= NMIN)
    k = p[p.pasa]
    n1, n2 = int(k["2e10"].sum()), int(k["1e12"].sum())
    if n1 and n2:
        r = n1 / n2; e = r * math.sqrt(1 / n1 + 1 / n2)
    else:
        r = e = float("nan")
    return p, n1, n2, r, e, list(k.index)


tab = pd.read_csv("tabla_regimen_pool.csv", dtype={"nota": str}).fillna({"nota": ""})
bak = pd.read_csv("tabla_regimen_pool_threshold001.csv")
targets = [("Zr", "liquido"), ("Ag", "liquido"), ("Cu", "liquido"), ("Cu", "vitreo")]
changes = []
for sp, reg in targets:
    p, n1, n2, r, e, kept = pooled(sp, reg)
    print(f"\n=== {sp} {reg} (ico default) -- T {WIN[reg][0]}..{WIN[reg][-1]}, conteos 2e10/1e12 por T ===")
    print("  " + "  ".join(f"{t}:{row['2e10']}/{row['1e12']}{'' if row.pasa else '×'}" for t, row in p.iterrows()))
    dropped = [t for t in WIN[reg] if t not in kept]
    print(f"  T descartados (N<{NMIN}): {dropped if dropped else 'ninguno'} | retenidos {len(kept)}/{len(WIN[reg])}"
          f" | pool {n1}/{n2} | razon {r:.4f} +- {e:.4f}")
    o = bak[(bak.especie == sp) & (bak.metodo == "Voronoi") & (bak.regimen == reg)].iloc[0]
    changes.append((sp, reg, o.N_pool_2e10, o.N_pool_1e12, o.razon, o.error, n1, n2, r, e, kept, dropped))
    if WRITE:
        m = (tab.especie == sp) & (tab.metodo == "Voronoi") & (tab.regimen == reg)
        assert m.sum() == 1
        i = tab.index[m][0]
        for c in ["N_pool_2e10", "N_pool_1e12", "razon", "error"]:
            tab[c] = tab[c].astype(object)
        if kept:
            for c, v in zip(["N_pool_2e10", "N_pool_1e12", "razon", "error"], [n1, n2, round(r, 4), round(e, 4)]):
                tab.at[i, c] = v
            tab.at[i, "nota"] = (f"ico default; T usados {' '.join(map(str, kept))}"
                                 + (f"; descartados por N<{NMIN}: {' '.join(map(str, dropped))}" if dropped else "")
                                 + f"; N>={NMIN} por T")
        else:
            for c in ["N_pool_2e10", "N_pool_1e12", "razon", "error"]:
                tab.at[i, c] = ""
            tab.at[i, "nota"] = "NO CALCULABLE: N<%d en todos los T (2e10/1e12: %s)" % (
                NMIN, "; ".join(f"T={t}: {row['2e10']}/{row['1e12']}" for t, row in p.iterrows()))

print("\n=== CAMBIOS vs version threshold=0.01 (incl. Zr/Ag vitreo de la Tarea O) ===")
print(f"{'fila':16s} {'antes (pool, razon+-err)':38s} {'ahora':38s}")
allrows = [("Zr", "vitreo"), ("Ag", "vitreo")] + targets
for sp, reg in allrows:
    o = bak[(bak.especie == sp) & (bak.metodo == "Voronoi") & (bak.regimen == reg)].iloc[0]
    p, n1, n2, r, e, kept = pooled(sp, reg)
    now = f"{n1}/{n2}  {r:.4f}+-{e:.4f}" if kept else "no calculable"
    print(f"{sp+'-'+reg:16s} {o.N_pool_2e10}/{o.N_pool_1e12}  {o.razon:.4f}+-{o.error:.4f}".ljust(56) + now)

if WRITE:
    tab.to_csv("tabla_regimen_pool.csv", index=False)
    print("\n-> tabla_regimen_pool.csv actualizada")
