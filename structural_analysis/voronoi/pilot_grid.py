import warnings; warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import numpy as np, sys
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier
RADII = {1: 1.60, 2: 1.28, 3: 1.44}; SP = {1: "Zr", 2: "Cu", 3: "Ag"}
def run(dump, radii, edge, rel):
    p = import_file(dump)
    if radii:
        pt = p.source.data.particles_.particle_types_
        for k, r in RADII.items(): pt.type_by_id_(k).radius = r
    p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True, use_radii=radii, edge_threshold=edge, relative_face_threshold=rel))
    d = p.compute(0)
    return d.particles["Coordination"][...], d.particles["Particle Type"][...], d.particles["Voronoi Index"][...]
dump = sys.argv[1]
for radii in (False, True):
    for edge, rel in ((0.0, 0.0), (0.1, 0.0), (0.0, 0.01), (0.1, 0.01)):
        c, t, vi = run(dump, radii, edge, rel)
        ico = ((vi[:,2]==0)&(vi[:,3]==0)&(vi[:,4]==12)&(vi[:,5]==0)).mean()*100
        s = "  ".join(f"{SP[k]} {c[t==k].mean():5.2f}(m{np.bincount(c[t==k]).argmax()})" for k in SP)
        print(f"radii={int(radii)} edge={edge:<3} rel={rel:<4} | CN global {c.mean():5.2f} | {s} | ico {ico:.2f}%")
