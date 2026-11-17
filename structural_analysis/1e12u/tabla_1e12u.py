#!/usr/bin/env python3
"""Tablas de estructura local con la 1e12 UNIFICADA (17/11/2026). Solo lee; escribe en 1e12u/.

Entradas: ../aggregate/tabla_compare_rates_3tasas.csv (2e10 y 1e11, sin cambios),
          abundance_voronoi_1e12u.csv (voronoi_1e12u.py) y ../boo/1e12u/NN.boo.dump (criterio ico de Ryltsev).
Salidas : tabla_compare_rates_u.csv (mismo formato: rate,T,especie,tipo_cluster,fraccion,N; rate = 2e10, 1e11, 1e12u)
          tabla_regimen_pool_u.csv (razones 2e10/1e12u y 1e11/1e12u por régimen, mismo método que Etapa 3:
          vítreo 500/400/300 K, líquido 1900-900 K, se descarta una T si alguna tasa tiene N < 20,
          razón = N_pool1/N_pool2, error = razón·sqrt(1/N1 + 1/N2)).
Chequeo : --validar reproduce con 2e10/1e12 (vieja) las razones de tabla_regimen_pool.csv y tabla_fk_zr.csv.
"""
import csv, sys
sys.path.insert(0, "../classify")
from boo_types import load_boo_dump, mask_ico, SPECIES  # noqa: E402

NESP = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
TIPOS = ["ico_voronoi", "ico_boo", "Z14_voronoi", "Z15_voronoi", "Z16_voronoi"]
VIT, LIQ = [500, 400, 300], list(range(1900, 800, -100))

tab = [r for r in csv.DictReader(open("../aggregate/tabla_compare_rates_3tasas.csv"))]


def N(rate, T, esp, tc, src):
    for r in src:
        if r["rate"] == rate and int(r["T"]) == T and r["especie"] == esp and r["tipo_cluster"] == tc:
            return int(r["N"])
    raise KeyError((rate, T, esp, tc))


def pool(a, b, esp, tc, Ts, src):
    kept = [T for T in Ts if N(a, T, esp, tc, src) >= 20 and N(b, T, esp, tc, src) >= 20]
    if not kept:
        return None
    n1 = sum(N(a, T, esp, tc, src) for T in kept); n2 = sum(N(b, T, esp, tc, src) for T in kept)
    r = n1 / n2
    return n1, n2, r, r * (1 / n1 + 1 / n2) ** 0.5, kept


if "--validar" in sys.argv:
    print("validación con 2e10/1e12 (vieja) contra tabla_regimen_pool.csv / tabla_fk_zr.csv:")
    for esp, tc, reg, Ts in (("Cu", "ico_boo", "vitreo", VIT), ("Cu", "ico_voronoi", "vitreo", VIT), ("Cu", "ico_voronoi", "liquido", LIQ),
                             ("Zr", "Z14_voronoi", "vitreo", VIT), ("Zr", "Z15_voronoi", "vitreo", VIT), ("Zr", "Z16_voronoi", "vitreo", VIT),
                             ("Ag", "ico_voronoi", "vitreo", VIT)):
        p = pool("2e10", "1e12", esp, tc, Ts, tab)
        print(f"  {esp} {tc:12s} {reg:8s} pool {p[0]}/{p[1]} razón {p[2]:.4f} ± {p[3]:.4f}")
    sys.exit()

# filas nuevas 1e12u
ab = {(int(r["T"]), r["especie"], r["tipo_cluster"]): float(r["fraccion"]) for r in csv.DictReader(open("abundance_voronoi_1e12u.csv"))}
nuevas = []
for nn in range(1, 18):
    T = 1900 - 100 * (nn - 1)
    for esp in ("Ag", "Cu", "Zr"):
        for tc in ("ico_voronoi", "Z14_voronoi", "Z15_voronoi", "Z16_voronoi"):
            f = ab[(T, esp, tc)]
            nuevas.append(dict(rate="1e12u", T=T, especie=esp, tipo_cluster=tc, fraccion=f"{f:.6f}", N=round(f * NESP[esp])))
    d = load_boo_dump(f"../boo/1e12u/{nn:02d}.boo.dump")
    types = d[:, 1].astype(int); ico = mask_ico(d)
    for t, esp in SPECIES.items():
        n = int((ico & (types == t)).sum())
        nuevas.append(dict(rate="1e12u", T=T, especie=esp, tipo_cluster="ico_boo", fraccion=f"{n / NESP[esp]:.6f}", N=n))
filas = [r for r in tab if r["rate"] in ("2e10", "1e11")] + nuevas
with open("tabla_compare_rates_u.csv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=["rate", "T", "especie", "tipo_cluster", "fraccion", "N"]); w.writeheader(); w.writerows(filas)

out = []
for a in ("2e10", "1e11"):
    for esp, tc in (("Cu", "ico_boo"), ("Cu", "ico_voronoi"), ("Zr", "Z14_voronoi"), ("Zr", "Z15_voronoi"), ("Zr", "Z16_voronoi"), ("Ag", "ico_voronoi")):   # Zr-ico no está en la tabla de Etapa 3 (N<20)
        for reg, Ts in (("vitreo", VIT), ("liquido", LIQ)):
            p = pool(a, "1e12u", esp, tc, Ts, filas)
            if p is None:
                out.append([f"{a}/1e12u", esp, tc, reg, "", "", "", "", "ninguna T con N>=20"]); continue
            n1, n2, r, e, kept = p
            out.append([f"{a}/1e12u", esp, tc, reg, n1, n2, f"{r:.4f}", f"{e:.4f}", " ".join(map(str, kept))])
            print(f"{a}/1e12u {esp} {tc:12s} {reg:8s} pool {n1:5d}/{n2:5d}  razón {r:.3f} ± {e:.3f}  ({len(kept)} T)")
with open("tabla_regimen_pool_u.csv", "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["comparacion", "especie", "tipo_cluster", "regimen", "N_pool_a", "N_pool_b", "razon", "error", "T_usados"]); w.writerows(out)
