#!/usr/bin/env python3
"""Tarea G - piloto de clasificacion Voronoi (2e10, config a 300K),
via OVITO (no LAMMPS nativo -- decision ya tomada por el usuario).

Lee snaps/17.dump de 2e10 (id type x y z, generado por
../boo/../../<rate>/snaps via restart_to_dump.lmp -- NO se toca ese
archivo ni structural_analysis/boo/ ni classify/, solo se lee).

Radios metalicos (Slater 1964, J. Chem. Phys. 41, 3199 -- confirmados
con el usuario 2026-09-17 antes de usarlos; Zr es HCP y Cu/Ag son FCC
en su forma pura, los tres nativamente CN=12, asi que el radio empirico
de Slater coincide exactamente con el radio de Goldschmidt CN12, sin
correccion de coordinacion):
  Zr = 1.60 A   Cu = 1.28 A   Ag = 1.44 A

Mapeo de tipos LAMMPS (pair_coeff * * ZrCuAg.lammps.eam Zr Cu Ag, mismo
usado en toda la Etapa 2): 1=Zr, 2=Cu, 3=Ag.

edge_threshold=0.1: es un valor ABSOLUTO en A (documentacion de OVITO:
"edges shorter than this threshold are ignored when counting the
number of edges of a face"), no una fraccion. Para este sistema la
distancia tipica primeros-vecinos es ~2.6-3.2 A (N=65536 en una caja de
~104.75 A de lado -> densidad ~0.0566 atomos/A^3 -> separacion tipica
(1/dens)^(1/3) ~ 2.65 A), asi que 0.1 A es ~4% de esa escala: solo
recorta aristas degeneradas de ruido numerico de la teselacion, no
aristas reales. Razonable para esta escala, se usa el valor de ejemplo
de la documentacion de OVITO tal cual.
"""
import numpy as np
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier

DUMP = "../2e10/snaps/17.dump"
RADII = {1: 1.60, 2: 1.28, 3: 1.44}   # Zr, Cu, Ag (Slater 1964)
SPECIES = {1: "Zr", 2: "Cu", 3: "Ag"}
EDGE_THRESHOLD = 0.1

pipeline = import_file(DUMP)
data = pipeline.compute(0)

print(f"=== {DUMP} ===")
print(f"frame timestep: {data.attributes.get('Timestep')}")
print(f"natoms: {data.particles.count}")
cell = data.cell
print(f"box (diag): {np.diag(cell[:, :3])}")

ptypes = data.particles_.particle_types_
for tid, name in SPECIES.items():
    ptypes.type_by_id_(tid).radius = RADII[tid]
    print(f"  set radius type {tid} ({name}) = {RADII[tid]} A")

modifier = VoronoiAnalysisModifier(
    compute_indices=True,
    use_radii=True,
    edge_threshold=EDGE_THRESHOLD,
)
pipeline.modifiers.append(modifier)
data = pipeline.compute(0)

print("\nparticle properties after Voronoi:", list(data.particles.keys()))

types = data.particles["Particle Type"][...]
coord = data.particles["Coordination"][...]
vindex = data.particles["Voronoi Index"][...]
vol = data.particles["Atomic Volume"][...]

print(f"\nVoronoi Index shape: {vindex.shape}  (componente i = # caras con i aristas)")
print(f"Coordination: min={coord.min()} max={coord.max()} mean={coord.mean():.3f}")
print(f"Atomic Volume: mean={vol.mean():.3f} A^3  (V_total/N = "
      f"{np.prod(np.diag(cell[:, :3]))/data.particles.count:.3f} A^3, deben coincidir)")

np.save("vindex_2e10_300K.npy", vindex)
np.save("types_2e10_300K.npy", types)
np.save("coord_2e10_300K.npy", coord)
np.save("vol_2e10_300K.npy", vol)
print("\n-> guardado vindex_2e10_300K.npy, types_2e10_300K.npy, "
      "coord_2e10_300K.npy, vol_2e10_300K.npy (para el resto de la Tarea G)")
