#!/usr/bin/env python3
"""Voronoi (OVITO VoronoiAnalysisModifier) sobre las 51 configs
(3 tasas x 17 T), snaps/NN.dump -> CN por especie + indices exactos.

Uso:  python voronoi_51.py [--radii]
  default : parametros default de OVITO (sin radios, sin umbrales de cara).
            Validado en piloto 2e10/300K contra Fig. 6 de Mendez & Jaroszewicz
            (CN global 13.99, Zr 15.4 / Cu 12.6 / Ag 14.0).
  --radii : sensibilidad, use_radii=True con radios de Slater (Zr 1.60, Cu 1.28,
            Ag 1.44 A) fijados en pipeline.source (NO sobre compute(0), que es
            una copia y no tiene efecto). Salidas con sufijo _radii.

Mapeo type->especie (pair_coeff * * ZrCuAg.lammps.eam Zr Cu Ag): 1=Zr 2=Cu 3=Ag.
Columnas de 'Voronoi Index' en OVITO 3.15: col i = caras con (i+1) aristas,
=> <n3,n4,n5,n6> = columnas 2..5 (las 0,1 son siempre 0).

NOTA Z16: el indice <0,0,12,4> es DEGENERADO en sistemas Cu-Zr (Lu et al.,
Phys. Rev. B 102, 184201, 2020): agrupa dos poliedros geometricamente
distintos y solo uno es el Z16 real. El Z16 por Voronoi es por tanto un
LIMITE SUPERIOR, a cruzar despues contra el Z16 de BOO.
"""
import sys, time, csv, warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import numpy as np
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier

USE_RADII = "--radii" in sys.argv
SUF = "_radii" if USE_RADII else ""
RATES = ["2e10", "1e11", "1e12"]
NCONF = 17
SP = {1: "Zr", 2: "Cu", 3: "Ag"}
RADII = {1: 1.60, 2: 1.28, 3: 1.44}
INDICES = {"ico_voronoi": (0, 0, 12, 0), "Z14_voronoi": (0, 0, 12, 2),
           "Z15_voronoi": (0, 0, 12, 3), "Z16_voronoi": (0, 0, 12, 4)}
NATOMS = 65536
CN_RANGE = range(6, 24)


def t_nominal(nn):
    return 1900 - (nn - 1) * 100


def voronoi_one(path):
    p = import_file(path)
    src = p.source.data
    assert src.particles.count == NATOMS, path
    assert tuple(src.cell.pbc) == (True, True, True), path
    if USE_RADII:
        pt = p.source.data.particles_.particle_types_
        for k, r in RADII.items():
            pt.type_by_id_(k).radius = r
    p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True, use_radii=USE_RADII))
    d = p.compute(0)
    vi = d.particles["Voronoi Index"][...]
    assert vi.shape[1] >= 6 and not vi[:, :2].any(), "convencion de columnas inesperada"
    return d.particles["Coordination"][...], d.particles["Particle Type"][...], vi[:, 2:6]


t_start = time.time()
cn_rows, ab_rows, summary = [], [], []
for rate in RATES:
    for nn in range(1, NCONF + 1):
        c, t, nv = voronoi_one(f"../{rate}/snaps/{nn:02d}.dump")
        T = t_nominal(nn)
        assert set(np.unique(t)) == set(SP), "tipos inesperados"
        groups = [("all", np.ones(len(t), bool))] + [(SP[k], t == k) for k in SP]
        for lab, m in groups:
            hist = np.bincount(c[m], minlength=max(CN_RANGE) + 1)
            assert hist[max(CN_RANGE) + 1:].sum() == 0 and hist[:min(CN_RANGE)].sum() == 0
            for n in CN_RANGE:
                cn_rows.append((rate, T, lab, n, hist[n] / m.sum()))
        for name, idx in INDICES.items():
            hit = (nv == np.array(idx)).all(axis=1)
            for k, sp in SP.items():
                m = t == k
                ab_rows.append((rate, T, sp, name, (hit & m).sum() / m.sum()))
        summary.append((rate, T, c.mean(), *[c[t == k].mean() for k in SP]))
    print(f"[{rate}] done  ({time.time()-t_start:.0f} s)", flush=True)

with open(f"cn_distribution_voronoi{SUF}.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "T", "especie", "CN", "fraccion"])
    w.writerows((r, T, s, n, f"{x:.6f}") for r, T, s, n, x in cn_rows)
with open(f"abundance_voronoi{SUF}.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "T", "especie", "tipo_cluster", "fraccion"])
    w.writerows((r, T, s, tc, f"{x:.6f}") for r, T, s, tc, x in ab_rows)

print(f"\nmodo: {'use_radii Slater' if USE_RADII else 'default OVITO'}  -- total {time.time()-t_start:.0f} s")
print("\nCN medio (global / Zr / Cu / Ag) por tasa y T:")
for r, T, *m in summary:
    print(f"  {r:>5} {T:5d}  " + "  ".join(f"{x:6.2f}" for x in m))
print("\nAVISO: <0,0,12,4> es degenerado en Cu-Zr (Lu et al., PRB 102, 184201, 2020); "
      "Z16_voronoi es un LIMITE SUPERIOR a cruzar contra Z16 de BOO.")
