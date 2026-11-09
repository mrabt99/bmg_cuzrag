#!/usr/bin/env python3
"""Voronoi (OVITO VoronoiAnalysisModifier, parametros DEFAULT: sin radios, sin umbrales)
sobre los snaps de qgrad_ext (qgrad_ext/snaps/<label>.dump). COPIA de qgrad/aggregate/voronoi_qgrad.py (que no se toca). Copia de la logica de
structural_analysis/voronoi/voronoi_51.py, que no se toca.

Uso (desde qgrad_ext/aggregate/): python voronoi_qgrad.py <label> [<label> ...]
Salidas (append-safe: se reescriben por completo con las labels dadas):
  abundance_voronoi_ext.csv       rate,T,especie,tipo_cluster,fraccion   (mismo formato que abundance_voronoi.csv)
  cn_distribution_voronoi_ext.csv rate,T,especie,CN,fraccion
  voronoi_atoms/<label>.npz         por atomo: id, type, coord, vindex(n3..n6)   (para la busqueda de Laves)
El T nominal viene de la etiqueta: NN=06..16 -> 1900-100*(NN-1); 17u/17r -> 300 K.
"""
import sys, csv, warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import numpy as np
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier

SP = {1: "Zr", 2: "Cu", 3: "Ag"}
INDICES = {"ico_voronoi": (0, 0, 12, 0), "Z14_voronoi": (0, 0, 12, 2),
           "Z15_voronoi": (0, 0, 12, 3), "Z16_voronoi": (0, 0, 12, 4)}
NATOMS = 65536
CN_RANGE = range(6, 24)
SNAPS = "../snaps"


def t_nominal(label):
    return 300 if label.startswith("17") else 1900 - (int(label) - 1) * 100


def voronoi_one(path):
    p = import_file(path)
    assert p.source.data.particles.count == NATOMS, path
    assert tuple(p.source.data.cell.pbc) == (True, True, True), path
    p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True))   # DEFAULT
    d = p.compute(0)
    vi = d.particles["Voronoi Index"][...]
    assert vi.shape[1] >= 6 and not vi[:, :2].any(), "convencion de columnas inesperada"
    return (d.particles["Particle Identifier"][...], d.particles["Particle Type"][...],
            d.particles["Coordination"][...], vi[:, 2:6])


labels = sys.argv[1:]
cn_rows, ab_rows = [], []
for lab in labels:
    ids, t, c, nv = voronoi_one(f"{SNAPS}/{lab}.dump")
    o = np.argsort(ids); ids, t, c, nv = ids[o], t[o], c[o], nv[o]
    T = t_nominal(lab)
    assert set(np.unique(t)) == set(SP)
    for lb, m in [("all", np.ones(len(t), bool))] + [(SP[k], t == k) for k in SP]:
        hist = np.bincount(c[m], minlength=max(CN_RANGE) + 1)
        assert hist[max(CN_RANGE) + 1:].sum() == 0 and hist[:min(CN_RANGE)].sum() == 0
        for n in CN_RANGE:
            cn_rows.append(("qgrad_ext", lab, T, lb, n, hist[n] / m.sum()))
    for name, idx in INDICES.items():
        hit = (nv == np.array(idx)).all(axis=1)
        for k, sp in SP.items():
            m = t == k
            ab_rows.append(("qgrad_ext", lab, T, sp, name, (hit & m).sum() / m.sum()))
    np.savez_compressed(f"voronoi_atoms/{lab}.npz", id=ids, type=t, coord=c, vindex=nv)
    print(f"{lab} T={T}: CN medio {c.mean():.2f}  ico_voro={int((nv == np.array((0,0,12,0))).all(1).sum())}", flush=True)

with open("abundance_voronoi_ext.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "config", "T", "especie", "tipo_cluster", "fraccion"])
    w.writerows((r, l, T, s, tc, f"{x:.6f}") for r, l, T, s, tc, x in ab_rows)
with open("cn_distribution_voronoi_ext.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "config", "T", "especie", "CN", "fraccion"])
    w.writerows((r, l, T, s, n, f"{x:.6f}") for r, l, T, s, n, x in cn_rows)
