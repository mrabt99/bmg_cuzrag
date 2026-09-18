#!/usr/bin/env python3
"""Tarea H, paso 2 -- Voronoi (ico estricto + Z14/Z15/Z16 por CN) sobre
los 51 dumps de posiciones (structural_analysis/<tasa>/snaps/NN.dump),
SOLO LECTURA de snaps/, no se toca boo/ ni voronoi/ (el piloto).

Decision de disenio final (Tarea H, ya no la de la Tarea F):
  n_ico       -> Voronoi <0,0,12,0>, columnas (2,3,4,5) del vector
                 'Voronoi Index' == (0,0,12,0). Convencion de columnas
                 confirmada empiricamente en el piloto (Tarea G):
                 columna i = caras con (i+1) aristas.
  Z14/Z15/Z16 -> Coordination == 14 / 15 / 16 (NO percentil de q4,
                 ese criterio de la Tarea F queda descartado).
Parametros de la teselacion: use_radii=True con Zr=1.60/Cu=1.28/Ag=1.44
A (Slater 1964, confirmados con el usuario en la Tarea G),
edge_threshold=0.1 A, relative_face_threshold=0.01 (nombre de keyword
confirmado por introspeccion de OVITO 3.15.5 instalado, Tarea G) --
ya verificado en el piloto que este filtro apenas cambia el resultado
(2992->2988 icosaedros en el caso de prueba), se deja puesto por ser el
criterio final acordado.
"""
import csv
import numpy as np
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier

RATES = ["2e10", "1e11", "1e12"]
NCONFIGS = 17
RADII = {1: 1.60, 2: 1.28, 3: 1.44}
SPECIES = {1: "Zr", 2: "Cu", 3: "Ag"}
EDGE_THRESHOLD = 0.1
RELATIVE_FACE_THRESHOLD = 0.01
N3, N4, N5, N6 = 2, 3, 4, 5


def t_nominal(nn):
    return 1900 - (nn - 1) * 100


def voronoi_one(dump_path):
    pipeline = import_file(dump_path)
    ptypes = pipeline.compute(0).particles_.particle_types_
    for tid in SPECIES:
        ptypes.type_by_id_(tid).radius = RADII[tid]
    pipeline.modifiers.append(VoronoiAnalysisModifier(
        compute_indices=True, use_radii=True,
        edge_threshold=EDGE_THRESHOLD,
        relative_face_threshold=RELATIVE_FACE_THRESHOLD,
    ))
    data = pipeline.compute(0)
    ids = data.particles["Particle Identifier"][...]
    order = np.argsort(ids)
    vindex = data.particles["Voronoi Index"][...][order]
    coord = data.particles["Coordination"][...][order]
    types = data.particles["Particle Type"][...][order]
    return vindex, coord, types


rows = []
for rate in RATES:
    for nn in range(1, NCONFIGS + 1):
        path = f"../{rate}/snaps/{nn:02d}.dump"
        vindex, coord, types = voronoi_one(path)

        strict_ico = ((vindex[:, N3] == 0) & (vindex[:, N4] == 0) &
                      (vindex[:, N5] == 12) & (vindex[:, N6] == 0))
        z14 = coord == 14
        z15 = coord == 15
        z16 = coord == 16

        T = t_nominal(nn)
        for t, name in SPECIES.items():
            m = types == t
            rows.append({
                "tasa": rate, "config": f"{nn:02d}", "T_nominal": T,
                "especie": name,
                "N_ico_Voronoi": int((strict_ico & m).sum()),
                "N_Z14": int((z14 & m).sum()),
                "N_Z15": int((z15 & m).sum()),
                "N_Z16": int((z16 & m).sum()),
            })
        print(f"[{rate}] {nn:02d} (T={T}K): ico_voro={strict_ico.sum()} "
              f"Z14={z14.sum()} Z15={z15.sum()} Z16={z16.sum()}")

with open("step2_voronoi.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"\n-> step2_voronoi.csv ({len(rows)} filas)")
