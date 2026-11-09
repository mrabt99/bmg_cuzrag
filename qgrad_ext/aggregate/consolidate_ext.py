#!/usr/bin/env python3
"""Tabla de clusters de la rama qgrad_ext (5 curvas: 2e10, 1e11, 1e12, qgrad, qgrad_ext). SOLO LECTURA de structural_analysis/;
escribe en qgrad_ext/aggregate/. Mismos criterios que qgrad/aggregate/consolidate_qgrad.py (BOO ico q6>=0.6 & W6hat<=-0.16 via classify/boo_types;
Voronoi OVITO DEFAULT, indices exactos <0,0,12,k>).
Curva qgrad_ext: configs 01-05 (1900-1500 K) = 2e10; 06-13 (1400-700 K) = qgrad (misma trayectoria, restart de 700 K identico);
14 (600), 15 (500), 16 (400), 17r (300 relajado), 17u (300 sin relajar) = NUEVAS, calculadas aca desde snaps/boo/voronoi_atoms de qgrad_ext.
Salida: tabla_compare_rates_5curvas.csv  rate,T,especie,tipo_cluster,fraccion,N,origen,config   (origen de la fila de qgrad_ext: 2e10|qgrad|ext)
        tabla_clusters_ext.csv            formato ancho solo para las configs nuevas
Uso (desde qgrad_ext/aggregate/): PYTHONDONTWRITEBYTECODE=1 python consolidate_ext.py"""
import sys, csv
SA = "/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis"
sys.path.insert(0, f"{SA}/classify")
import numpy as np
from boo_types import load_boo_dump, mask_ico

NSP = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
TYPE_OF = {"Zr": 1, "Cu": 2, "Ag": 3}
VORO_IDX = {"ico_voronoi": (0, 0, 12, 0), "Z14_voronoi": (0, 0, 12, 2), "Z15_voronoi": (0, 0, 12, 3), "Z16_voronoi": (0, 0, 12, 4)}
NEW = [("14", 600), ("15", 500), ("16", 400), ("17u", 300), ("17r", 300)]

rows = list(csv.DictReader(open(f"{SA}/qgrad/aggregate/tabla_compare_rates_4curvas.csv")))
out = [(r["rate"], int(r["T"]), r["especie"], r["tipo_cluster"], float(r["fraccion"]), int(r["N"]), "", r["config_qgrad"]) for r in rows]
# qgrad_ext = copia de las filas de qgrad con config 01..13
for r in rows:
    if r["rate"] == "qgrad" and r["config_qgrad"] and r["config_qgrad"][:2].isdigit() and 1 <= int(r["config_qgrad"][:2]) <= 13 and not r["config_qgrad"].startswith("17"):
        nn = int(r["config_qgrad"]); org = "2e10" if nn <= 5 else "qgrad"
        out.append(("qgrad_ext", int(r["T"]), r["especie"], r["tipo_cluster"], float(r["fraccion"]), int(r["N"]), org, r["config_qgrad"]))
wide = []
for lab, T in NEW:
    d = load_boo_dump(f"../boo/{lab}.boo.dump"); t = d[:, 1].astype(int); ico = mask_ico(d)
    z = np.load(f"voronoi_atoms/{lab}.npz"); tv, nv = z["type"], z["vindex"]
    assert len(t) == len(tv) == 65536 and (np.sort(z["id"]) == z["id"]).all()
    bc = {sp: int((ico & (t == TYPE_OF[sp])).sum()) for sp in NSP}
    vc = {(sp, tc): int((((nv == np.array(ix)).all(1)) & (tv == TYPE_OF[sp])).sum()) for sp in NSP for tc, ix in VORO_IDX.items()}
    for sp in ("Ag", "Cu", "Zr"):
        wide.append(("qgrad_ext", lab, T, sp, NSP[sp], bc[sp], vc[(sp, "ico_voronoi")], vc[(sp, "Z14_voronoi")], vc[(sp, "Z15_voronoi")], vc[(sp, "Z16_voronoi")]))
    out.append(("qgrad_ext", T, "Cu", "ico_boo", bc["Cu"] / NSP["Cu"], bc["Cu"], "ext", lab))
    for (sp, tc), n in vc.items(): out.append(("qgrad_ext", T, sp, tc, n / NSP[sp], n, "ext", lab))
    print(f"{lab} T={T}: Cu-ico BOO={bc['Cu']}  Cu-ico Voro={vc[('Cu','ico_voronoi')]}  Zr Z14/Z15/Z16={vc[('Zr','Z14_voronoi')]}/{vc[('Zr','Z15_voronoi')]}/{vc[('Zr','Z16_voronoi')]}")
ORDER = ["2e10", "1e11", "1e12", "qgrad", "qgrad_ext"]
out.sort(key=lambda r: (r[2], r[3], ORDER.index(r[0]), -r[1], r[7]))
with open("tabla_compare_rates_5curvas.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "T", "especie", "tipo_cluster", "fraccion", "N", "origen", "config"])
    w.writerows((a, b, c, d, f"{e:.6f}", n, o, l) for a, b, c, d, e, n, o, l in out)
with open("tabla_clusters_ext.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["tasa", "config", "T_nominal", "especie", "N_total_especie", "N_ico_BOO", "N_ico_Voronoi", "N_Z14", "N_Z15", "N_Z16"]); w.writerows(wide)
print(f"tabla_compare_rates_5curvas.csv: {len(out)} filas")
