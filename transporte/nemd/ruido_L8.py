"""Ruido real de kappa NEMD en L = 419 A, Q = 0.5 eV/ps (tag L8, tramos 1-6 = 6 ns; se descarta 1 ns).
Pregunta: el error por superbloques de 500 ps (+-0.007 con 2 ns) estaba subestimado? El ruido estimado por el modo termico
mas lento daba +-0.017 para 2 ns (-> ~+-0.011 para 5 ns si escala con 1/sqrt(t)).
Mismo procedimiento que analiza_nemd.py (perfil medio, recta en el 60 % central de cada mitad, J = Q/(2A)).
1) verificacion de la serie (bloques de 100 ps contiguos, sin repetidos)
2) kappa con 1-3 ns (como antes), 3-6 ns (nuevo) y 1-6 ns
3) analisis de bloqueo: sem con superbloques de m = 1, 2, 5, 10, 25 bloques (100 ps - 2.5 ns)
4) autocorrelacion lag-1 de kappa por bloque de 100 ps y por superbloque de 500 ps
Uso (desde nemd/): python ruido_L8.py -> ruido_L8.txt, ruido_L8.png"""
import glob, re, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Q=0.5; W=10.0; VMEAN=143891.93343795; A=VMEAN**(2/3); CONV=1.602176634e-19/1e-12/1e-10/2; TEQ=1000.0; FR=0.6
GK=(0.779,0.017); EST2NS=0.017   # GK mismo vidrio; ruido estimado (modo lento) para 2 ns
def leer(tag):
    blk=[]
    for f in sorted(glob.glob(f"{tag}_p*_Tprof.dat"),key=lambda s:int(re.search(r"_p(\d+)_",s).group(1))):
        L=[l.split() for l in open(f) if not l.startswith("#")]; i=0
        while i<len(L):
            ts,nc=int(L[i][0]),int(L[i][1]); d=np.array(L[i+1:i+1+nc],float); blk.append((ts,d[:,1],d[:,3])); i+=1+nc
    box=[l for l in open(sorted(glob.glob(f"log.{tag}_p1.*"))[0]) if "orthogonal box" in l][-1]
    nums=[float(x) for x in re.findall(r"[-0-9.eE+]+",box.split("=")[1])]; return blk,nums[2],nums[5]-nums[2]
def kappa(z,T,z0,Lz):
    s=(z-(z0+W/2))%Lz; h=Lz/2; a=(1-FR)/2*h; g=[]
    for lo in (0,h):
        m=(s>lo+a)&(s<lo+h-a); g.append(abs(np.polyfit(s[m],T[m],1)[0]))
    return CONV*Q/(A*np.mean(g))
def ac1(x):
    x=np.asarray(x)-np.mean(x); return np.sum(x[1:]*x[:-1])/np.sum(x*x)
out=[]
blk,z0,Lz=leer("L8"); ts=np.array([b[0] for b in blk]); dts=np.diff(ts)
out.append(f"L8: {len(blk)} bloques, pasos {ts[0]}..{ts[-1]}, saltos entre bloques: {sorted(set(dts.tolist()))} (esperado [100000])")
t=ts*1e-3; z=blk[0][1]; P=np.array([b[2] for b in blk]); sel=t>TEQ; Ps=P[sel]; ns=len(Ps)
kb=np.array([kappa(z,p,z0,Lz) for p in Ps])   # kappa por bloque de 100 ps
def tramo(lo,hi):
    m=(t>lo)&(t<=hi); Pm=P[m]; k=kappa(z,Pm.mean(0),z0,Lz); nsb=len(Pm)//5
    ksb=[kappa(z,Pm[i*5:(i+1)*5].mean(0),z0,Lz) for i in range(nsb)]
    return k,np.std(ksb,ddof=1)/np.sqrt(nsb),nsb
for lo,hi,lab in ((1000,3000,"1-3 ns (datos de antes)"),(3000,6000,"3-6 ns (extension)"),(1000,6000,"1-6 ns (todo)")):
    k,e,nsb=tramo(lo,hi); out.append(f"{lab:26s} kappa = {k:.3f} +- {e:.3f} (sem de {nsb} superbloques de 500 ps)")
kA,eA,_=tramo(1000,3000); kB,eB,_=tramo(3000,6000)
out.append(f"diferencia (3-6) - (1-3) ns = {kB-kA:+.3f}  (con sem de superbloques: +- {np.hypot(eA,eB):.3f};"
           f" con ruido estimado {EST2NS:.3f}*sqrt(2/t): +- {np.hypot(EST2NS,EST2NS*np.sqrt(2/3)):.3f})")
ktot=kappa(z,Ps.mean(0),z0,Lz)
out.append(""); out.append(f"analisis de bloqueo ({ns} bloques de 100 ps, kappa del perfil medio = {ktot:.3f}):")
out.append("  m  tam.(ps)  n_sb   sem      sem/kappa   lag-1 entre superbloques")
MS=[1,2,5,10,25]; sems=[]
for m in MS:
    nsb=ns//m; ksb=np.array([kappa(z,Ps[i*m:(i+1)*m].mean(0),z0,Lz) for i in range(nsb)])
    e=np.std(ksb,ddof=1)/np.sqrt(nsb); sems.append((m,e,e/np.sqrt(2*(nsb-1)))); r=ac1(ksb) if nsb>3 else np.nan
    out.append(f"{m:3d} {m*100:8d} {nsb:5d}   {e:.4f}   {100*e/ktot:5.1f} %     {r:+.2f}")
out.append(f"lag-1 de kappa por bloque de 100 ps: {ac1(kb):+.2f}   (Q = 0.25 daba +0.59)")
out.append(f"ruido estimado (modo lento) escalado a {ns*0.1:.0f} ns: {EST2NS*np.sqrt(2/(ns*0.1)):.4f}")
out.append(f"referencia GK, mismo vidrio: {GK[0]} +- {GK[1]}")
fig,ax=plt.subplots(1,2,figsize=(10,3.8))
ax[0].plot(t[sel],kb,".-",lw=0.6,label="bloques de 100 ps")
nsb=ns//5; ax[0].plot(t[sel][4::5][:nsb],[kappa(z,Ps[i*5:(i+1)*5].mean(0),z0,Lz) for i in range(nsb)],"o-",label="superbloques de 500 ps")
ax[0].axhline(ktot,color="C3",lw=0.8,label=f"1-6 ns: {ktot:.3f}"); ax[0].axhline(GK[0],color="k",lw=0.8,label="GK")
ax[0].axvline(3000,color="gray",ls=":",lw=0.8); ax[0].set_xlabel("t (ps)"); ax[0].set_ylabel(r"$\kappa$ (W/m/K)"); ax[0].legend(fontsize=7)
m_,e_,ee_=map(np.array,zip(*sems)); ax[1].errorbar(m_*100,e_,ee_,fmt="o-")
ax[1].axhline(EST2NS*np.sqrt(2/(ns*0.1)),color="C3",ls="--",lw=0.8,label="ruido estimado (modo lento)")
ax[1].set_xscale("log"); ax[1].set_xlabel("tamano de superbloque (ps)"); ax[1].set_ylabel(r"sem de $\kappa$ (W/m/K)"); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig("ruido_L8.png",dpi=150)
open("ruido_L8.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
