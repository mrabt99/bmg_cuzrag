"""Analisis NEMD (fix ehex) de kappa a 300 K. Uso (desde nemd/): python analiza_nemd.py [n ...]   (default 4 8 16)
Lee L<n>_p*_Tprof.dat (bloques de 100 ps), descarta el transitorio TEQ[n] ps (0 si n no esta en TEQ), promedia el perfil T(z) y ajusta
rectas en la parte central (FR = fraccion de cada mitad) de las dos mitades fria->caliente y caliente->fria.
kappa = J/|dT/dz|, J = Q/(2A); en unidades: kappa[W/m/K] = 801.088*Q[eV/ps]/(A[A^2]*g[K/A]).
Error: dispersion entre super-bloques de 500 ps del regimen estacionario. Sistematico: rango de ajuste 0.5/0.6/0.7.
Extrapolacion 1/kappa = 1/kappa_inf + c/L si hay >= 2 longitudes. Salidas: nemd.txt, nemd.png"""
import sys, glob, re, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Q=0.5; W=10.0; VMEAN=143891.93343795; A=VMEAN**(2/3); CONV=1.602176634e-19/1e-12/1e-10/2   # J=Q/(2A), g en K/A
TEQ={4:500.0,8:1000.0,16:3000.0}; FRS=(0.5,0.6,0.7); SB=5   # super-bloque = 5 bloques de 100 ps
NS=[int(x) for x in sys.argv[1:]] or [4,8,16]
GK=(0.779,0.017)   # cal_300K_1e12_N8192, mismo vidrio (trapecio, 0.5-2 ps)
def leer(n):
    blk=[]
    for f in sorted(glob.glob(f"L{n}_p*_Tprof.dat"),key=lambda s:int(re.search(r"_p(\d+)_",s).group(1))):
        L=[l.split() for l in open(f) if not l.startswith("#")]; i=0
        while i<len(L):
            ts,nc=int(L[i][0]),int(L[i][1]); d=np.array(L[i+1:i+1+nc],float); blk.append((ts,d[:,1],d[:,3])); i+=1+nc
    log=sorted(glob.glob(f"log.L{n}_p1.*"))[0]
    box=[l for l in open(log) if "orthogonal box" in l][-1]
    nums=[float(x) for x in re.findall(r"[-0-9.eE+]+",box.split("=")[1])]
    return blk,nums[2],nums[5]-nums[2]
def kappa(z,T,z0,Lz,fr):
    s=(z-(z0+W/2))%Lz; h=Lz/2; a=(1-fr)/2*h; g=[]
    for lo in (0,h):
        m=(s>lo+a)&(s<lo+h-a); g.append(abs(np.polyfit(s[m],T[m],1)[0]))
    return CONV*Q/(A*np.mean(g)), g
out=[]; res={}
fig,ax=plt.subplots(1,3,figsize=(14,3.8))
for n in NS:
    blk,z0,Lz=leer(n); t=np.array([b[0] for b in blk])*1e-3   # ps
    kb=np.array([kappa(b[1],b[2],z0,Lz,0.6)[0] for b in blk])
    ax[1].plot(t,kb,".-",label=f"L={Lz:.0f} A")
    sel=[b for b,tt in zip(blk,t) if tt>TEQ.get(n,0.0)]
    if len(sel)<SB: out.append(f"n={n}: solo {len(sel)} bloques estacionarios, insuficiente"); continue
    z=sel[0][1]; Tm=np.mean([b[2] for b in sel],0)
    ks=[kappa(z,Tm,z0,Lz,fr)[0] for fr in FRS]; k=ks[1]; g=kappa(z,Tm,z0,Lz,0.6)[1]
    nsb=len(sel)//SB; ksb=[kappa(z,np.mean([b[2] for b in sel[i*SB:(i+1)*SB]],0),z0,Lz,0.6)[0] for i in range(nsb)]
    e=np.std(ksb,ddof=1)/np.sqrt(nsb) if nsb>1 else np.nan
    res[n]=(Lz,k,e,np.std(ks,ddof=1))
    out.append(f"n={n:2d} L={Lz:6.1f} A  bloques est.={len(sel):3d}  kappa={k:.3f} +- {e:.3f} (est, {nsb} superbloques)"
               f" +- {np.std(ks,ddof=1):.3f} (rango ajuste {FRS}: {', '.join(f'{x:.3f}' for x in ks)})"
               f"  |dT/dz| mitades = {g[0]:.4f}, {g[1]:.4f} K/A  Tmin/Tmax = {Tm.min():.1f}/{Tm.max():.1f} K")
    ax[0].plot((z-z0)/Lz,Tm,label=f"L={Lz:.0f} A")
if len(res)>=2:
    x=np.array([1/v[0] for v in res.values()]); y=np.array([1/v[1] for v in res.values()])
    wts=np.array([v[1]**2/max(v[2],1e-3) for v in res.values()])   # 1/sigma de 1/kappa
    p=np.polyfit(x,y,1,w=wts); out.append(f"extrapolacion 1/kappa vs 1/L: kappa_inf = {1/p[1]:.3f} W/m/K (pendiente {p[0]:.2f} A m K/W)")
    ax[2].errorbar(x,y,[v[2]/v[1]**2 for v in res.values()],fmt="o"); xx=np.linspace(0,x.max()*1.1,10); ax[2].plot(xx,np.polyval(p,xx),"--")
    ax[2].axhline(1/GK[0],color="k",lw=0.8,label="GK (8192 at.)")
out.append(f"referencia GK, mismo vidrio: {GK[0]} +- {GK[1]} W/m/K")
ax[0].set_xlabel("z/L"); ax[0].set_ylabel("T (K)"); ax[0].legend(fontsize=8)
ax[1].axhline(GK[0],color="k",lw=0.8); ax[1].set_xlabel("t (ps)"); ax[1].set_ylabel(r"$\kappa$ por bloque de 100 ps"); ax[1].legend(fontsize=8)
ax[2].set_xlabel("1/L (1/A)"); ax[2].set_ylabel(r"1/$\kappa$ (m K/W)"); ax[2].set_xlim(left=0); ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig("nemd.png",dpi=150)
open("nemd.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
