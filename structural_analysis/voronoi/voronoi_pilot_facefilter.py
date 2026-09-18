#!/usr/bin/env python3
"""Tarea G, chequeo puntual -- mismo piloto que voronoi_pilot.py (2e10,
config a 300K), agregando filtro de area de cara.

Nombre exacto del keyword confirmado por introspeccion de la version
instalada (OVITO 3.15.5, `help(VoronoiAnalysisModifier)`), no asumido
de la documentacion web: `relative_face_threshold` (float, fraccion del
area TOTAL de la superficie del poliedro; 0.01 = descarta caras con
menos del 1% del area total). Se aplica simultaneamente con
`face_threshold` (area absoluta, default 0 = sin efecto, no se toca).

edge_threshold se deja igual (0.1 A, ya justificado en voronoi_pilot.py).
"""
import numpy as np
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier

DUMP = "../2e10/snaps/17.dump"
RADII = {1: 1.60, 2: 1.28, 3: 1.44}
SPECIES = {1: "Zr", 2: "Cu", 3: "Ag"}
EDGE_THRESHOLD = 0.1
RELATIVE_FACE_THRESHOLD = 0.01

pipeline = import_file(DUMP)
ptypes = pipeline.compute(0).particles_.particle_types_
for tid in SPECIES:
    ptypes.type_by_id_(tid).radius = RADII[tid]

modifier = VoronoiAnalysisModifier(
    compute_indices=True,
    use_radii=True,
    edge_threshold=EDGE_THRESHOLD,
    relative_face_threshold=RELATIVE_FACE_THRESHOLD,
)
pipeline.modifiers.append(modifier)
data = pipeline.compute(0)

print(f"=== {DUMP}  (relative_face_threshold={RELATIVE_FACE_THRESHOLD}, "
      f"edge_threshold={EDGE_THRESHOLD}) ===")
coord = data.particles["Coordination"][...]
vindex = data.particles["Voronoi Index"][...]
ids = data.particles["Particle Identifier"][...]
types = data.particles["Particle Type"][...]

order = np.argsort(ids)
np.save("vindex_2e10_300K_facefilter.npy", vindex[order])
np.save("types_2e10_300K_facefilter.npy", types[order])
np.save("coord_2e10_300K_facefilter.npy", coord[order])
np.save("ids_2e10_300K_facefilter.npy", ids[order])

print(f"Coordination: min={coord.min()} max={coord.max()} mean={coord.mean():.3f}")
print(f"vindex shape: {vindex.shape}")
print("-> vindex/types/coord/ids _facefilter.npy guardados")
