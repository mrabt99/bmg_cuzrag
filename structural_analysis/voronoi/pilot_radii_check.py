import warnings; warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import numpy as np
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier, ComputePropertyModifier
RADII = {1: 1.60, 2: 1.28, 3: 1.44}; SP = {1: "Zr", 2: "Cu", 3: "Ag"}
DUMP = "../2e10/snaps/17.dump"

def vor(p):
    p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True, use_radii=True, edge_threshold=0.1, relative_face_threshold=0.01))
    d = p.compute(0); return d

def report(lab, d):
    c = d.particles["Coordination"][...]; t = d.particles["Particle Type"][...]; vi = d.particles["Voronoi Index"][...]
    ico = ((vi[:,2]==0)&(vi[:,3]==0)&(vi[:,4]==12)&(vi[:,5]==0)).sum()
    print(f"\n[{lab}] CN medio {c.mean():.3f}  ico {ico}")
    print("    CN " + "".join(f"{n:6d}" for n in range(9,19)))
    for k in SP:
        m = t == k
        print(f"{SP[k]:>6} " + "".join(f"{(c[m]==n).mean():6.3f}" for n in range(9,19)) + f"  mean={c[m].mean():.2f} mode={np.bincount(c[m]).argmax()}")

# A: metodo del script anterior (radios sobre copia de compute) -> se espera que NO tenga efecto
p = import_file(DUMP); pt = p.compute(0).particles_.particle_types_
for k, r in RADII.items(): pt.type_by_id_(k).radius = r
report("A: radios sobre copia de compute(0)", vor(p))

# B: radios en la fuente del pipeline
p = import_file(DUMP)
pt = p.source.data.particles_.particle_types_
for k, r in RADII.items(): pt.type_by_id_(k).radius = r
report("B: radios en pipeline.source", vor(p))

# C: propiedad Radius por particula via ComputePropertyModifier (inequivoco)
p = import_file(DUMP)
p.modifiers.append(ComputePropertyModifier(output_property="Radius",
    expressions=["ParticleType==1 ? 1.60 : (ParticleType==2 ? 1.28 : 1.44)"]))
d = vor(p)
r = d.particles["Radius"][...]; t = d.particles["Particle Type"][...]
print("\nRadius por tipo (C):", {SP[k]: float(np.unique(r[t==k])[0]) for k in SP})
report("C: Radius por particula", d)
