"""Prueba de tamaño: kappa con 8192 vs 65536 átomos, mismo criterio que mail_seba/kappa_vs_T.py.
Trapecio, bloques de 100 ps (10 por réplica de 1 ns), ventana [0.5,2] ps con test de pendiente
(|pendiente| < 2*error bootstrap), sistemático = std de ventanas [0.5,1],[1,2],[0.5,1.5],[0.75,2].
Además kappa(t) en 5, 10, 20 ps (cola larga, hipótesis de propagones/tamaño).
Uso (desde transporte/): python size_8192/compara_tamano.py -> size_8192/compara_tamano.{txt,png}"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
kB=1.380649e-23; eV=1.602176634e-19; NB=10; TMAX=20.0
W0,W1=0.5,2.0; WINS=[(0.5,1),(1,2),(0.5,1.5),(0.75,2)]
CASOS=[("cal_300K_1e12",300,3,"N=65536, 1e12"),("cal_300K_1e12_N8192",300,3,"N=8192, 1e12"),
       ("cal_1000K_2e10",1000,6,"N=65536, 2e10 (6 rép.)"),("cal_1000K_1e12_N8192",1000,3,"N=8192, 1e12")]
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
def blocks(d,T,nrep):
    V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); K=[]
    for r in range(1,nrep+1):
        J=np.loadtxt(f"{d}/cal_r{r}_hflux.dat",comments="#"); dt=(J[1,0]-J[0,0])*1e-3; nl=int(TMAX/dt)+1
        X=[J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)]; m=len(X[0])//NB
        B=np.array([np.mean([acf(x[b*m:(b+1)*m],nl) for x in X],axis=0) for b in range(NB)])
        K.append((np.cumsum(B,1)-0.5*B[:,:1])*dt/(V*1e-30*kB*T**2)*1e-12)
    return np.concatenate(K),dt
rng=np.random.default_rng(0); out=[]; res={}
out.append(f"{'caso':24s} {'kappa':>6s} {'sem':>6s} {'sist':>6s} {'pend':>8s} {'+-boot':>7s} meseta  k(5ps)        k(10ps)       k(20ps)")
fig,axs=plt.subplots(1,2,figsize=(10,3.8))
for d,T,nrep,lab in CASOS:
    K,dt=blocks(d,T,nrep); t=np.arange(K.shape[1])*dt; sel=(t>=W0)&(t<=W1); nb=len(K)
    g=K[:,sel].mean(1); k=g.mean(); e=g.std(ddof=1)/np.sqrt(nb)
    sist=np.std([K[:,(t>=a)&(t<=b)].mean() for a,b in WINS],ddof=1)
    sl=lambda M: np.polyfit(t[sel],M[sel],1)[0]
    s0=sl(K.mean(0)); se=np.std([sl(K[rng.integers(0,nb,nb)].mean(0)) for _ in range(400)])
    lt=[(K[:,int(round(x/dt))].mean(),K[:,int(round(x/dt))].std(ddof=1)/np.sqrt(nb)) for x in (5,10,20)]
    res[d]=(k,e,sist)
    out.append(f"{lab:24s} {k:6.3f} {e:6.3f} {sist:6.3f} {s0:+8.4f} {se:7.4f} {'si' if abs(s0)<2*se else 'NO':6s} "
               +"  ".join(f"{a:.3f}±{b:.3f}" for a,b in lt))
    ax=axs[0 if T==300 else 1]; m=K.mean(0); s=K.std(0,ddof=1)/np.sqrt(nb)
    ax.plot(t,m,label=lab); ax.fill_between(t,m-s,m+s,alpha=0.25)
for i,(a,b) in enumerate([("cal_300K_1e12","cal_300K_1e12_N8192"),("cal_1000K_2e10","cal_1000K_1e12_N8192")]):
    da=res[b][0]-res[a][0]; ea=np.hypot(res[a][1],res[b][1])
    out.append(f"diferencia {b} - {a}: {da:+.3f} ± {ea:.3f} ({abs(da)/ea:.1f} sigma, solo error estadístico)")
for ax,T in zip(axs,(300,1000)):
    ax.set_xscale("log"); ax.set_xlim(0.02,TMAX); ax.axvspan(W0,W1,color="0.9"); ax.set_title(f"{T} K")
    ax.set_xlabel("t (ps)"); ax.set_ylabel(r"$\kappa(t)$ (W/m/K)"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("size_8192/compara_tamano.png",dpi=150)
open("size_8192/compara_tamano.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
