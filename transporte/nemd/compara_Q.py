"""Prueba de linealidad NEMD: kappa a L = 419 A con Q = 0.5 (tag L8) y Q = 0.25 eV/ps (tag L8q025).
Mismo procedimiento que analiza_nemd.py: descarta 1 ns, perfil medio, recta en el 60 % central de cada mitad,
error por superbloques de 500 ps, sistematico por rango de ajuste (0.5/0.6/0.7).
Ademas: curvatura del perfil (ajuste cuadratico en cada mitad) como indicador directo de no linealidad.
Uso (desde nemd/): python compara_Q.py -> compara_Q.txt, compara_Q.png"""
import glob, re, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
W=10.0; VMEAN=143891.93343795; A=VMEAN**(2/3); CONV=1.602176634e-19/1e-12/1e-10/2; TEQ=1000.0; SB=5; FRS=(0.5,0.6,0.7)
CASOS=[("L8",0.5),("L8q025",0.25)]
def leer(tag):
    blk=[]
    for f in sorted(glob.glob(f"{tag}_p*_Tprof.dat"),key=lambda s:int(re.search(r"_p(\d+)_",s).group(1))):
        L=[l.split() for l in open(f) if not l.startswith("#")]; i=0
        while i<len(L):
            ts,nc=int(L[i][0]),int(L[i][1]); d=np.array(L[i+1:i+1+nc],float); blk.append((ts,d[:,1],d[:,3])); i+=1+nc
    box=[l for l in open(sorted(glob.glob(f"log.{tag}_p1.*"))[0]) if "orthogonal box" in l][-1]
    nums=[float(x) for x in re.findall(r"[-0-9.eE+]+",box.split("=")[1])]; return blk,nums[2],nums[5]-nums[2]
def ajuste(z,T,z0,Lz,fr,Q,deg=1):
    s=(z-(z0+W/2))%Lz; h=Lz/2; a=(1-fr)/2*h; g=[]; c=[]
    for lo in (0,h):
        m=(s>lo+a)&(s<lo+h-a); p=np.polyfit(s[m]-(lo+h/2),T[m],deg); g.append(abs(p[-2])); c.append(p[0] if deg==2 else 0)
    return CONV*Q/(A*np.mean(g)), g, c
out=[]; res={}; fig,ax=plt.subplots(1,2,figsize=(10,3.8))
for tag,Q in CASOS:
    blk,z0,Lz=leer(tag); t=np.array([b[0] for b in blk])*1e-3; sel=[b for b,tt in zip(blk,t) if tt>TEQ]
    z=sel[0][1]; Tm=np.mean([b[2] for b in sel],0)
    ks=[ajuste(z,Tm,z0,Lz,fr,Q)[0] for fr in FRS]; k=ks[1]
    nsb=len(sel)//SB; ksb=[ajuste(z,np.mean([b[2] for b in sel[i*SB:(i+1)*SB]],0),z0,Lz,0.6,Q)[0] for i in range(nsb)]
    e=np.std(ksb,ddof=1)/np.sqrt(nsb); sist=np.std(ks,ddof=1)
    _,g2,c2=ajuste(z,Tm,z0,Lz,0.6,Q,deg=2)
    res[tag]=(k,e,sist)
    out.append(f"{tag:7s} Q={Q:.2f} eV/ps  bloques est.={len(sel):3d}  kappa = {k:.3f} +- {e:.3f} (est, {nsb} superbloques) +- {sist:.3f} (ajuste)"
               f"  dT/dz = {np.mean(ajuste(z,Tm,z0,Lz,0.6,Q)[1]):.4f} K/A  curvatura d2T/dz2 = {c2[0]*2:+.2e}, {c2[1]*2:+.2e} K/A^2"
               f"  Tmin/Tmax = {Tm.min():.1f}/{Tm.max():.1f}")
    ax[0].plot((z-z0)/Lz,Tm-Tm.mean(),label=f"Q={Q} eV/ps"); ax[0].plot((z-z0)/Lz,(Tm-Tm.mean())*0.5/Q,"--",lw=0.8,label=f"Q={Q} escalado a 0.5")
    kb=[ajuste(b[1],b[2],z0,Lz,0.6,Q)[0] for b in blk]; ax[1].plot(t,kb,".-",label=f"Q={Q}")
(k1,e1,s1),(k2,e2,s2)=res["L8"],res["L8q025"]; d=k2-k1; ed=np.hypot(e1,e2)
out.append(f"diferencia kappa(Q=0.25) - kappa(Q=0.5) = {d:+.3f} +- {ed:.3f} ({abs(d)/ed:.1f} sigma, estadistico)")
ax[0].set_xlabel("z/L"); ax[0].set_ylabel("T - <T> (K)"); ax[0].legend(fontsize=7)
ax[1].axhline(0.779,color="k",lw=0.8); ax[1].set_ylim(0,2); ax[1].set_xlabel("t (ps)"); ax[1].set_ylabel(r"$\kappa$ por bloque de 100 ps"); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig("compara_Q.png",dpi=150)
open("compara_Q.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
