import warnings; warnings.filterwarnings("ignore")
import numpy as np, csv
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier
SP={1:"Zr",2:"Cu",3:"Ag"}
rows=[]; bad=0
for rate in ["2e10","1e11","1e12"]:
    for nn in range(1,18):
        p=import_file(f"../../{rate}/snaps/{nn:02d}.dump")
        p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True))
        d=p.compute(0)
        c=d.particles["Coordination"][...]; t=d.particles["Particle Type"][...]; vi=d.particles["Voronoi Index"][...]
        T=1900-(nn-1)*100
        base=(vi[:,2]==0)&(vi[:,3]==0)&(vi[:,4]==12)
        big=vi[:,6:].sum(axis=1)  # caras con >=7 aristas
        for k,s in SP.items():
            m=t==k; n=m.sum()
            for j,name in [(0,"ico"),(2,"Z14"),(3,"Z15"),(4,"Z16")]:
                e51=(base&(vi[:,5]==j)&m); e_strict=(e51&(big==0))
                bad+=int(e51.sum()!=e_strict.sum())
                cn=(c==(12+j))&m if j else (c==12)&m
                # CN==12+j SIN indice: poblacion laxa; interseccion con indice exacto
                rows.append((rate,T,s,name,e51.sum()/n,e_strict.sum()/n,cn.sum()/n,(cn&e51).sum()/max(cn.sum(),1)))
print("celdas donde n>=7 != 0 cambia el conteo exacto:",bad)
with open("a2_recalculo_51configs.csv","w") as f:
    w=csv.writer(f); w.writerow(["rate","T","especie","tc","exact_cols25","exact_strict","cn_default","frac_cn_que_es_exacto"]); w.writerows(rows)
