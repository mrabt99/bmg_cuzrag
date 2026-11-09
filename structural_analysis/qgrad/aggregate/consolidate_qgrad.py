#!/usr/bin/env python3
"""Etapa 3 (paso 1) para qgrad. Lee SOLO LECTURA las tablas/dumps de las 3 tasas originales y las de qgrad, escribe en qgrad/aggregate/:
  tabla_clusters_qgrad.csv           formato ancho como tabla_clusters_51.csv (tasa,config,T_nominal,especie,N_total_especie,
                                     N_ico_BOO,N_ico_Voronoi,N_Z14,N_Z15,N_Z16) + compartido. Voronoi = OVITO DEFAULT
                                     (tabla_clusters_51.csv usa use_radii+umbral de cara; NO son comparables entre si).
  tabla_compare_rates_4curvas.csv    formato tidy como tabla_compare_rates_3tasas.csv (rate,T,especie,tipo_cluster,fraccion,N)
                                     + compartido (1 = fila de qgrad identica a 2e10, T>=1500 K) + config (qgrad).
Configs de qgrad: 01-05 (1900-1500 K) = compartidas con 2e10 (misma trayectoria; no se recalculan);
06-16 = seg01..seg11; 300 K = 17r (relajado, si existe) y/o 17u (sin relajar).
Uso: python consolidate_qgrad.py   (desde qgrad/aggregate/; PYTHONDONTWRITEBYTECODE=1)"""
import sys, os, csv
sys.path.insert(0, "../../classify")
import numpy as np
from boo_types import load_boo_dump, mask_ico, SPECIES

NSP = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
TYPE_OF = {"Zr": 1, "Cu": 2, "Ag": 3}
VORO_IDX = {"ico_voronoi": (0, 0, 12, 0), "Z14_voronoi": (0, 0, 12, 2), "Z15_voronoi": (0, 0, 12, 3), "Z16_voronoi": (0, 0, 12, 4)}
ORDER = ["2e10", "1e11", "1e12", "qgrad"]
T_OF = lambda nn: 1900 - (nn - 1) * 100

# --- fuentes de las 3 tasas (solo lectura)
ab = {}
for r in csv.DictReader(open("../../voronoi/abundance_voronoi.csv")):
    ab[(r["rate"], int(r["T"]), r["especie"], r["tipo_cluster"])] = round(float(r["fraccion"]) * NSP[r["especie"]])
boo_cu = {}
for r in csv.DictReader(open("../../aggregate/tabla_compare_rates_3tasas.csv")):
    if r["tipo_cluster"] == "ico_boo":
        boo_cu[(r["rate"], int(r["T"]))] = int(r["N"])

def boo_counts(path):
    d = load_boo_dump(path); t = d[:, 1].astype(int); ico = mask_ico(d)
    return {sp: int((ico & (t == TYPE_OF[sp])).sum()) for sp in NSP}

# --- configs de qgrad
cfgs = []   # (label, T, compartido, boo_path, voronoi_getter)
for nn in range(1, 6):
    cfgs.append((f"{nn:02d}", T_OF(nn), 1, f"../../boo/2e10/{nn:02d}.boo.dump", ("2e10", T_OF(nn))))
for nn in range(6, 17):
    cfgs.append((f"{nn:02d}", T_OF(nn), 0, f"../boo/{nn:02d}.boo.dump", f"voronoi_atoms/{nn:02d}.npz"))
for lab in ("17u", "17r"):
    if os.path.exists(f"../boo/{lab}.boo.dump") and os.path.exists(f"voronoi_atoms/{lab}.npz"):
        cfgs.append((lab, 300, 0, f"../boo/{lab}.boo.dump", f"voronoi_atoms/{lab}.npz"))

def voro_counts(src):
    out = {}
    if isinstance(src, tuple):                 # compartido: se toma de la tabla default de 2e10
        rate, T = src
        for sp in NSP:
            for tc in VORO_IDX: out[(sp, tc)] = ab[(rate, T, sp, tc)]
        return out
    z = np.load(src); t, nv = z["type"], z["vindex"]
    for sp in NSP:
        for tc, idx in VORO_IDX.items():
            out[(sp, tc)] = int((((nv == np.array(idx)).all(1)) & (t == TYPE_OF[sp])).sum())
    return out

wide, tidy_q = [], []
for lab, T, comp, bpath, vsrc in cfgs:
    bc = boo_counts(bpath); vc = voro_counts(vsrc)
    for sp in ("Ag", "Cu", "Zr"):
        wide.append(("qgrad", lab, T, sp, NSP[sp], bc[sp], vc[(sp, "ico_voronoi")], vc[(sp, "Z14_voronoi")],
                     vc[(sp, "Z15_voronoi")], vc[(sp, "Z16_voronoi")], comp))
    tidy_q.append((lab, T, comp, "Cu", "ico_boo", bc["Cu"]))
    for (sp, tc), n in vc.items(): tidy_q.append((lab, T, comp, sp, tc, n))
# chequeo: los compartidos (2e10) reproducen el BOO ya tabulado
for lab, T, comp, bpath, _ in cfgs:
    if comp: assert boo_counts(bpath)["Cu"] == boo_cu[("2e10", T)], f"BOO compartido {lab} no coincide con la tabla 3 tasas"

with open("tabla_clusters_qgrad.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["tasa", "config", "T_nominal", "especie", "N_total_especie", "N_ico_BOO", "N_ico_Voronoi", "N_Z14", "N_Z15", "N_Z16", "compartido_con_2e10"])
    w.writerows(wide)

rows = []   # rate,T,especie,tipo_cluster,fraccion,N,compartido,config
for r in csv.DictReader(open("../../aggregate/tabla_compare_rates_3tasas.csv")):
    rows.append((r["rate"], int(r["T"]), r["especie"], r["tipo_cluster"], float(r["fraccion"]), int(r["N"]), 0, ""))
have = {(r[0], r[1], r[2], r[3]) for r in rows}
for (rate, T, sp, tc), n in ab.items():          # combos Voronoi que la tabla 3 tasas no incluia (Zr-ico, Ag-FK, Cu-FK...)
    if (rate, T, sp, tc) not in have: rows.append((rate, T, sp, tc, n / NSP[sp], n, 0, ""))
for lab, T, comp, sp, tc, n in tidy_q:
    variante = "qgrad"
    rows.append((variante, T, sp, tc, n / NSP[sp], n, comp, lab))
rows.sort(key=lambda r: (r[2], r[3], ORDER.index(r[0]), -r[1], r[7]))
with open("tabla_compare_rates_4curvas.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "T", "especie", "tipo_cluster", "fraccion", "N", "compartido_con_2e10", "config_qgrad"])
    w.writerows((a, b, c, d, f"{e:.6f}", n, s, l) for a, b, c, d, e, n, s, l in rows)
print(f"tabla_clusters_qgrad.csv: {len(wide)} filas; tabla_compare_rates_4curvas.csv: {len(rows)} filas")
print("configs qgrad:", [(c[0], c[1], "compartido" if c[2] else "nuevo") for c in cfgs])
