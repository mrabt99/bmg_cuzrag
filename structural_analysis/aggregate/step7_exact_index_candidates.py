#!/usr/bin/env python3
"""Tarea N (pasos 1-4) -- Z14/Z15/Z16 por INDICE EXACTO <0,0,12,k> (OVITO default,
sin threshold), desde ../voronoi/abundance_voronoi.csv. NO usa CN== en ningun paso.

1. Adaptador: dtype rate=str, quita sufijo _voronoi de tipo_cluster, chequea 1:1.
2. Chequeo: Zr 300K/2e10 reproduce 1.6/2.8/1.5 %.
3. Conteos N = fraccion * N_total_especie (N_total constante: Zr 31129, Cu 31130,
   Ag 3277, de tabla_clusters_51.csv); filtro N>=20 en AMBAS tasas por punto de T
   antes de agregar; ventana vitrea COMPARTIDA (interseccion 2e10/1e12).
4. Muestra los 3 candidatos por especie (conteos, sin razones) y se detiene.

Nota (no bloqueante): Z16 exacto <0,0,12,4> es degenerado en Cu-Zr (Lu et al.,
Phys. Rev. B 102, 184201, 2020) => limite superior; ambiguedad no resuelta aqui.
"""
import pandas as pd

NTOT = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
NMIN = 20
RATES = ["2e10", "1e12"]
TG = {"2e10": (646.1 + 637.2) / 2, "1e12": (703.2 + 697.9) / 2}   # igual que step6
BUFFER = 100.0

# ---- paso 1: adaptador
raw = pd.read_csv("../voronoi/abundance_voronoi.csv", dtype={"rate": str})
tid = raw.copy()
tid["tipo_cluster"] = tid["tipo_cluster"].str.removesuffix("_voronoi")
assert len(raw) == len(tid) == 612
assert not tid.fraccion.isna().any()
assert not tid.duplicated(["rate", "T", "especie", "tipo_cluster"]).any(), "duplicados"
assert (tid.groupby(["rate", "T", "especie", "tipo_cluster"]).size() == 1).all()
print(f"paso 1: {len(raw)} filas -> {len(tid)} filas adaptadas, claves unicas, sin NaN")
print("  tipo_cluster:", sorted(raw.tipo_cluster.unique()), "->", sorted(tid.tipo_cluster.unique()))
print("  rate:", sorted(tid.rate.unique()), " T:", tid["T"].min(), "-", tid["T"].max(), f"({tid['T'].nunique()} T)")
tid.to_csv("abundance_voronoi_adaptado.csv", index=False)

# ---- conteos crudos (fraccion*Ntot debe ser entero)
tid["N"] = tid.fraccion * tid.especie.map(NTOT)
resid = (tid.N - tid.N.round()).abs().max()
bound = 0.5e-6 * max(NTOT.values())   # fraccion se guardo con 6 decimales
assert resid <= bound + 1e-9, f"fraccion*Ntot no entero (max resid {resid} > {bound})"
tid["N"] = tid.N.round().astype(int)
print(f"  conteos recuperados de fraccion*N_total: max |residuo a entero| = {resid:.4f}")

# ---- paso 2
chk = tid[(tid.rate == "2e10") & (tid["T"] == 300) & (tid.especie == "Zr")].set_index("tipo_cluster")
print("\npaso 2: Zr, 2e10, 300 K, indice exacto (crudo del CSV):")
for tc in ["Z14", "Z15", "Z16"]:
    print(f"  {tc}: fraccion={100*chk.loc[tc,'fraccion']:.3f} %  N={chk.loc[tc,'N']}")
# valores impresos en el piloto (2 decimales, 2026-09-18): 1.63/2.80/1.45.
# (el resumen escrito decia "1.5" para Z16: era 1.449 mal redondeado, no un dato distinto)
rep = {"Z14": 1.63, "Z15": 2.80, "Z16": 1.45}
ok = all(abs(100 * chk.loc[tc, "fraccion"] - v) <= 0.006 for tc, v in rep.items())
print("  impreso en el piloto 1.63/2.80/1.45 ->", "COINCIDE (a 2 decimales)" if ok else "NO COINCIDE")
assert ok, "no reproduce los numeros reportados: parar"

# ---- paso 3: ventana vitrea compartida y filtro N>=20
Ts = {r: sorted(tid[tid.rate == r]["T"].unique(), reverse=True) for r in RATES}
vit = {r: {t for t in Ts[r] if t < TG[r] - BUFFER} for r in RATES}
shared = sorted(vit["2e10"] & vit["1e12"], reverse=True)
print(f"\npaso 3: vitreo propio 2e10={sorted(vit['2e10'],reverse=True)} 1e12={sorted(vit['1e12'],reverse=True)}")
print(f"  ventana COMPARTIDA = {shared}")

piv = tid[tid["T"].isin(shared)].pivot_table(index=["especie", "tipo_cluster", "T"], columns="rate", values="N")
piv = piv[RATES].astype(int)
piv["pasa_N>=20"] = (piv["2e10"] >= NMIN) & (piv["1e12"] >= NMIN)

print(f"\npaso 4: candidatos por especie -- conteos crudos por T, ventana {shared} (filtro N>={NMIN} en ambas tasas)")
for sp in ["Zr", "Ag"]:
    print(f"\n=== {sp} (N_total {NTOT[sp]}) ===")
    print(f"{'cluster':8s} " + "".join(f"{'T='+str(t):>20s}" for t in shared) + f"{'N_pool 2e10':>14s}{'N_pool 1e12':>14s}{'T que pasan':>14s}")
    for tc in ["Z14", "Z15", "Z16"]:
        s = piv.loc[(sp, tc)]
        kept = s[s["pasa_N>=20"]]
        cells = "".join(f"{str(s.loc[t,'2e10'])+' / '+str(s.loc[t,'1e12'])+('' if s.loc[t,'pasa_N>=20'] else ' x'):>20s}" for t in shared)
        print(f"{tc:8s} {cells}{kept['2e10'].sum():14d}{kept['1e12'].sum():14d}{len(kept):>10d}/{len(shared)}")
print("\n(celda = N_2e10 / N_1e12; 'x' = descartado por N<20 en alguna tasa; pool = suma de T que pasan)")
