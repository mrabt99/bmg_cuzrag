import warnings; warnings.filterwarnings("ignore")
import numpy as np, collections
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier
SP={1:"Zr",2:"Cu",3:"Ag"}; RAD={1:1.60,2:1.28,3:1.44}
def run(path, mode):
    p=import_file(path)
    if mode=="A_bug":  # como step2: radios sobre compute(0) (copia) + umbrales
        pt=p.compute(0).particles_.particle_types_
        for k,r in RAD.items(): pt.type_by_id_(k).radius=r
        p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True,use_radii=True,edge_threshold=0.1,relative_face_threshold=0.01))
    elif mode=="thr_only":
        p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True,edge_threshold=0.1,relative_face_threshold=0.01))
    elif mode=="default":
        p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True))
    d=p.compute(0)
    return d.particles["Coordination"][...], d.particles["Particle Type"][...], d.particles["Voronoi Index"][...]
def stats(c,t,vi):
    out={}
    for k,s in SP.items():
        m=t==k; n=m.sum()
        ex={j:(((vi[:,2]==0)&(vi[:,3]==0)&(vi[:,4]==12)&(vi[:,5]==j)&(vi[:,6:].sum(axis=1)==0))&m).sum() for j in (0,2,3,4)}
        # nota: indice <0,0,12,k> exige n7+ = 0 tambien; comparar con el criterio de voronoi_51 (solo cols 2..5)
        ex51={j:(((vi[:,2]==0)&(vi[:,3]==0)&(vi[:,4]==12)&(vi[:,5]==j))&m).sum() for j in (0,2,3,4)}
        cn={j:((c==j)&m).sum() for j in (12,14,15,16)}
        out[s]=(n,ex,ex51,cn)
    return out
import sys
for rate,nn in [("2e10",17),("2e10",9),("2e10",1),("1e12",17),("1e11",17)]:
    path=f"../../{rate}/snaps/{nn:02d}.dump"
    print(f"\n##### {rate} config {nn:02d} (T={1900-(nn-1)*100}K)")
    for mode in ["default","thr_only","A_bug"]:
        c,t,vi=run(path,mode)
        st=stats(c,t,vi)
        print(f" [{mode}] CN medio={c.mean():.3f}  vi.shape={vi.shape}")
        for s in ("Zr","Cu","Ag"):
            n,ex,ex51,cn=st[s]
            print(f"   {s} n={n:5d} | exacto(n>=7=0) Z14/15/16={ex[2]}/{ex[3]}/{ex[4]}  exacto(cols2-5) {ex51[2]}/{ex51[3]}/{ex51[4]} ico={ex51[0]} | CN==14/15/16={cn[14]}/{cn[15]}/{cn[16]}")
