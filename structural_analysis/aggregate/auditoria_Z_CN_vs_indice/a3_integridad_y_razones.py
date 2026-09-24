import pandas as pd, numpy as np, math
NT={"Zr":31129,"Cu":31130,"Ag":3277}
# integridad de la cadena step2 -> tabla_clusters_51 -> abundancia
s2=pd.read_csv("../step2_voronoi.csv",dtype={"tasa":str,"config":str}); tb=pd.read_csv("../tabla_clusters_51.csv",dtype={"tasa":str,"config":str})
m=s2.merge(tb,on=["tasa","config","T_nominal","especie"],suffixes=("_s2","_tb"))
print("step2 vs tabla_clusters_51: filas",len(s2),len(tb),len(m),"; difs N_Z14/15/16:",[int((m[f"N_{z}_s2"]!=m[f"N_{z}_tb"]).sum()) for z in ("Z14","Z15","Z16")])
print("N_total_especie unicos:",tb.groupby("especie").N_total_especie.unique().to_dict())
bad=0
for r in ["2e10","1e11","1e12"]:
    a=pd.read_csv(f"../abundancia_{r}.csv")
    t=tb[tb.tasa==r].rename(columns={"T_nominal":"T"}).merge(a,on=["T","especie"])
    for z in ("Z14","Z15","Z16"):
        bad+=int(((t[f"N_{z}"]/t.N_total_especie)-t[f"n_{z}"]).abs().max()>5e-7)
print("abundancia vs tabla (N/Ntot) discrepancias:",bad)
# razones vitreas 2e10/1e12 (T=500,400,300; N>=20 por T) segun definicion
f=pd.read_csv("a2_recalculo_51configs.csv",dtype={"rate":str})
def ratio(get,sp,tc):
    p=[]
    for T in (500,400,300):
        n1=get("2e10",T,sp,tc); n2=get("1e12",T,sp,tc)
        if n1>=20 and n2>=20: p.append((n1,n2))
    if not p: return None
    n1=sum(a for a,_ in p); n2=sum(b for _,b in p); r=n1/n2
    return round(r,3),round(r*math.sqrt(1/n1+1/n2),3),n1,n2,len(p)
def g_ex(r,T,sp,tc): x=f[(f.rate==r)&(f["T"]==T)&(f.especie==sp)&(f.tc==tc)].exact_cols25.iloc[0]; return round(x*NT[sp])
def g_cnd(r,T,sp,tc): x=f[(f.rate==r)&(f["T"]==T)&(f.especie==sp)&(f.tc==tc)].cn_default.iloc[0]; return round(x*NT[sp])
def g_cnt(r,T,sp,tc):
    a=pd.read_csv(f"../abundancia_{r}.csv"); x=a[(a["T"]==T)&(a.especie==sp)][f"n_{tc}"].iloc[0]; return round(x*NT[sp])
print("\nRazon vitrea 2e10/1e12 (500/400/300K), razon ± error, N pool, T que pasan N>=20:")
for sp in ("Zr","Cu","Ag"):
    for tc in ("Z14","Z15","Z16"):
        print(f" {sp} {tc}: exacto={ratio(g_ex,sp,tc)} | CN default={ratio(g_cnd,sp,tc)} | CN+umbral(figura actual)={ratio(g_cnt,sp,tc)}")
