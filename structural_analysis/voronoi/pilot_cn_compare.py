#!/usr/bin/env python3
"""Piloto Voronoi (2e10, 300 K): default vs use_radii (Slater) -- CN global y
por especie, para comparar contra la Fig. 6 de Mendez & Jaroszewicz."""
import time, warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import numpy as np
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier

DUMP = "../2e10/snaps/17.dump"
RADII = {1: 1.60, 2: 1.28, 3: 1.44}   # Slater 1964
SP = {1: "Zr", 2: "Cu", 3: "Ag"}

def run(use_radii):
    t0 = time.time()
    p = import_file(DUMP)
    if use_radii:
        pt = p.compute(0).particles_.particle_types_
        for k, r in RADII.items():
            pt.type_by_id_(k).radius = r
    p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True, use_radii=use_radii,
                                               edge_threshold=0.1, relative_face_threshold=0.01))
    d = p.compute(0)
    dt = time.time() - t0
    return d.particles["Coordination"][...], d.particles["Particle Type"][...], d.particles["Voronoi Index"][...], dt

for ur in (False, True):
    c, t, vi, dt = run(ur)
    print(f"\n=== use_radii={ur}  ({dt:.1f} s) ===  CN medio global {c.mean():.3f}")
    rng = range(8, 20)
    print(f"{'CN':>6}" + "".join(f"{n:6d}" for n in rng))
    for lab, m in [("all", np.ones_like(t, bool))] + [(SP[k], t == k) for k in SP]:
        f = [(c[m] == n).mean() for n in rng]
        print(f"{lab:>6}" + "".join(f"{x:6.3f}" for x in f) + f"   mean={c[m].mean():.2f} mode={np.bincount(c[m]).argmax()}")
    ico = (vi[:,2]==0)&(vi[:,3]==0)&(vi[:,4]==12)&(vi[:,5]==0)
    print("ico <0,0,12,0>:", ico.sum(), f"({ico.mean()*100:.2f} %)")
