#!/usr/bin/env python3
"""Tarea O -- (1) filas Zr/Ag (vitreo, ventana compartida 500/400/300 K) de
tabla_regimen_pool.csv recalculadas con ico_voronoi DEFAULT (OVITO sin threshold,
../voronoi/abundance_voronoi.csv) en vez de la version con relative_face_threshold=0.01;
(2) FK exacto de Zr (<0,0,12,2/3/4>) en tabla_fk_zr.csv, aparte.

Reglas: rate=str; filtro N>=20 POR PUNTO DE T en ambas tasas (T que falla se
descarta del pool); razon agregada N_pool(2e10)/N_pool(1e12), error
razon*sqrt(1/N1+1/N2) sobre los conteos sumados. Sin CN== en ningun calculo.
Uso: python step8_ico_default_and_fk_zr.py [--write]
"""
import sys, math
import pandas as pd

WRITE = "--write" in sys.argv
NTOT = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
NMIN = 20
SHARED = [500, 400, 300]

a = pd.read_csv("../voronoi/abundance_voronoi.csv", dtype={"rate": str})
a["tipo_cluster"] = a.tipo_cluster.str.removesuffix("_voronoi")
a["N"] = (a.fraccion * a.especie.map(NTOT)).round().astype(int)
assert ((a.fraccion * a.especie.map(NTOT)) - a.N).abs().max() <= 0.5e-6 * 31130 + 1e-9


def counts(sp, tc):
    s = a[(a.especie == sp) & (a.tipo_cluster == tc) & a["T"].isin(SHARED) & a.rate.isin(["2e10", "1e12"])]
    p = s.pivot(index="T", columns="rate", values="N").loc[SHARED]
    p["pasa"] = (p["2e10"] >= NMIN) & (p["1e12"] >= NMIN)
    return p


def agg(p):
    k = p[p.pasa]
    n1, n2 = int(k["2e10"].sum()), int(k["1e12"].sum())
    if n1 == 0 or n2 == 0:
        return n1, n2, float("nan"), float("nan"), list(k.index)
    r = n1 / n2
    return n1, n2, r, r * math.sqrt(1 / n1 + 1 / n2), list(k.index)


# --- chequeo previo: ico default 300K/2e10 (Zr, Ag)
print("chequeo ico default, 2e10, 300 K (crudo del CSV):")
for sp in ["Zr", "Ag"]:
    r = a[(a.rate == "2e10") & (a["T"] == 300) & (a.especie == sp) & (a.tipo_cluster == "ico")].iloc[0]
    print(f"  {sp}: fraccion={100*r.fraccion:.3f} %  N={r.N}")

old = pd.read_csv("tabla_regimen_pool.csv")
newrows, old_vs_new = [], []
print(f"\n=== TAREA 1: ico default, conteos por T (N_2e10 / N_1e12), filtro N>={NMIN} por T ===")
for sp in ["Zr", "Ag"]:
    p = counts(sp, "ico")
    n1, n2, r, e, kept = agg(p)
    print(f"\n{sp}:")
    for t, row in p.iterrows():
        print(f"  T={t}: {row['2e10']:4d} / {row['1e12']:4d}  {'ok' if row.pasa else 'DESCARTADO (N<20)'}")
    print(f"  T retenidos: {kept}  pool {n1} / {n2}  razon {r:.4f} +- {e:.4f}" if kept else f"  NINGUN T pasa el filtro -> fila no calculable")
    o = old[(old.especie == sp) & (old.regimen == "vitreo")].iloc[0]
    print(f"  version con threshold=0.01: pool {o.N_pool_2e10} / {o.N_pool_1e12}  razon {o.razon:.4f} +- {o.error:.4f}")
    newrows.append(dict(especie=sp, metodo="Voronoi", regimen="vitreo", N_pool_2e10=n1, N_pool_1e12=n2,
                        razon=round(r, 4), error=round(e, 4), T_usados=" ".join(map(str, kept))))

print(f"\n=== TAREA 2: FK exacto Zr, conteos por T (N_2e10 / N_1e12) ===")
fk = []
for tc in ["Z14", "Z15", "Z16"]:
    p = counts("Zr", tc)
    n1, n2, r, e, kept = agg(p)
    print(f"  {tc}: " + "  ".join(f"T={t}: {row['2e10']}/{row['1e12']}{'' if row.pasa else ' x'}" for t, row in p.iterrows())
          + f"  | pool {n1}/{n2}  razon {r:.4f} +- {e:.4f}")
    fk.append(dict(especie="Zr", tipo_cluster=tc, regimen="vitreo", N_pool_2e10=n1, N_pool_1e12=n2,
                   razon=round(r, 4), error=round(e, 4), T_usados=" ".join(map(str, kept))))

if WRITE:
    import shutil, os
    bak = "tabla_regimen_pool_threshold001.csv"
    if not os.path.exists(bak):
        shutil.copy("tabla_regimen_pool.csv", bak)   # version previa (ico con threshold=0.01)
    tab = old.copy()
    tab["nota"] = ""
    for nr in newrows:
        m = (tab.especie == nr["especie"]) & (tab.metodo == "Voronoi") & (tab.regimen == "vitreo")
        assert m.sum() == 1
        i = tab.index[m][0]
        if nr["T_usados"] == "":
            p = counts(nr["especie"], "ico")
            det = "; ".join(f"T={t}: {r['2e10']}/{r['1e12']}" for t, r in p.iterrows())
            for c in ["N_pool_2e10", "N_pool_1e12", "razon", "error"]:
                tab[c] = tab[c].astype(object)
                tab.at[i, c] = ""
            tab.at[i, "nota"] = (f"NO CALCULABLE: ico default con N<{NMIN} en 2e10 en todos los T (N_2e10/N_1e12 -> {det}); "
                                 "filtro por T no forzado. Valor previo (threshold=0.01): 0.7159+-0.1182 descartado")
        else:
            for c in ["N_pool_2e10", "N_pool_1e12", "razon", "error"]:
                tab.at[i, c] = nr[c]
            tab.at[i, "nota"] = f"ico default (OVITO sin threshold); T usados {nr['T_usados']}; N>={NMIN} por T"
    tab.to_csv("tabla_regimen_pool.csv", index=False)
    f = pd.DataFrame(fk)
    f["nota"] = ""
    f.loc[f.tipo_cluster == "Z16", "nota"] = ("Z16 exacto <0,0,12,4> es limite superior: degenerado en Cu-Zr "
                                              "(Lu et al., Phys. Rev. B 102, 184201, 2020)")
    f.to_csv("tabla_fk_zr.csv", index=False)
    print("\n-> tabla_regimen_pool.csv actualizada (respaldo: %s), tabla_fk_zr.csv escrita" % bak)
