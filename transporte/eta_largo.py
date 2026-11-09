"""(1) eta a 800 y 900 K con las corridas largas (long_r4..r6, tramos de 1 ns concatenados, se descarta la fila repetida).
    Bloques de BL ps (500 ps), ACF de pxy,pxz,pyz hasta TMAX ps, trapecio. Se informa eta(t), ACF, meseta (pendiente
    con bootstrap en [t_a,t_b]) y regla tipo Zhang (std/media entre bloques = 30/40/50 %, valor en [tc/2,tc]).
    Tambien kappa (0.5-2 ps) y D de las corridas largas.
(2) kappa a 1000 K con 6 replicas (cal_r1..r6): ventana 0.5-2 ps y cola 2-20 ps.
Figura eta_largo.png. Uso: python eta_largo.py"""
import numpy as np, glob, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
kB=1.380649e-23; eV=1.602176634e-19; rng=np.random.default_rng(1)
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
def cat(d,tag,kind,npart):
    A=[np.loadtxt(f"{d}/{tag}_p{p}_{kind}.dat",comments="#") for p in range(1,npart+1)]
    return np.vstack([A[0]]+[a[1:] for a in A[1:]])
def blocks(cols,nl,m): 
    nb=len(cols[0])//m; return [np.mean([acf(c[b*m:(b+1)*m],nl) for c in cols],axis=0) for b in range(nb)]
sem=lambda a: a.std(ddof=1)/np.sqrt(len(a))
fig,ax=plt.subplots(1,3,figsize=(15,4.6)); a1,a2,a3=ax
CFG={800:(5,200.0,(50,150)),900:(3,100.0,(20,60))}; BL=500.0
for T,(npart,TMAX,(ta,tb)) in CFG.items():
    d=f"cal_{T}K_2e10"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); BP=[];BJ=[];D=[]
    for tag in ("long_r4","long_r5","long_r6"):
        P=cat(d,tag,"press",npart); dtp=(P[1,0]-P[0,0])*1e-3; m=int(BL/dtp)
        BP+=blocks([P[:,i]*1e5 for i in (4,5,6)],int(TMAX/dtp),m)
        J=cat(d,tag,"hflux",npart); dtj=(J[1,0]-J[0,0])*1e-3
        BJ+=blocks([J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)],int(3/dtj),int(BL/dtj))
        for p in range(1,npart+1):
            ms=np.loadtxt(f"{d}/{tag}_p{p}_msd.dat",comments="#"); t=ms[:,0]*1e-3; t=t-t[0]; h=len(t)//2
            D.append([np.polyfit(t[h:],ms[h:,i],1)[0]/6*10 for i in (1,2,3)])
    BP=np.array(BP); BJ=np.array(BJ); nb=len(BP)
    E=(np.cumsum(BP,1)-0.5*BP[:,:1])*dtp*V*1e-30/(kB*T)*1e-12*1e3; tp=np.arange(E.shape[1])*dtp
    K=(np.cumsum(BJ,1)-0.5*BJ[:,:1])*dtj/(V*1e-30*kB*T**2)*1e-12; tj=np.arange(K.shape[1])*dtj
    rho=BP.mean(0)/BP[:,0].mean()
    print(f"\n===== {T} K: {nb} bloques de {BL:.0f} ps (3 replicas x {npart} ns), ACF hasta {TMAX:.0f} ps")
    for x in (0.1,0.3,0.5,1,2): print(f"  ACF presion C/C0 a {x:>4} * TMAX ({x*TMAX:6.1f} ps): {rho[min(int(x*TMAX/dtp),len(rho)-1)]:+.3f}")
    print("  eta(t) [mPa s]: "+"  ".join(f"t={x:g}:{E[:,int(x/dtp)].mean():.1f}+-{sem(E[:,int(x/dtp)]):.1f}" for x in (5,10,20,50,100,150,199) if x<TMAX))
    sel=(tp>=ta)&(tp<=tb); g=E[:,sel].mean(1); sl=lambda M: np.polyfit(tp[sel],M[sel],1)[0]
    s0=sl(E.mean(0)); sb=np.std([sl(E[rng.integers(0,nb,nb)].mean(0)) for _ in range(300)])
    print(f"  ventana [{ta},{tb}] ps: eta = {g.mean():.1f} +- {sem(g):.1f} mPa s; pendiente {s0:+.3f} +- {sb:.3f} mPa s/ps -> {'meseta' if abs(s0)<2*sb else 'SIN meseta'}")
    r=E.std(0,ddof=1)/np.abs(E.mean(0)); out=[]
    for fr in (0.3,0.4,0.5):
        i=np.where((tp>=1)&(r>=fr))[0]; ic=int(i[0]) if len(i) else len(tp)-1; s=(tp>=tp[ic]/2)&(tp<=tp[ic])
        out.append(f"{100*fr:.0f}%: tc={tp[ic]:.0f} ps, eta={E[:,s].mean():.1f}")
    print("  regla Zhang: "+" | ".join(out))
    s5=(tj>=0.5)&(tj<=2); gk=K[:,s5].mean(1); Dm=np.mean(D,axis=0)
    print(f"  kappa[0.5-2 ps] = {gk.mean():.3f} +- {sem(gk):.3f} W/m/K ({nb} bloques) | D(Zr,Cu,Ag) = "+" ".join(f"{x:.3f}" for x in Dm)+" 1e-9 m2/s")
    e=E.mean(0); es=E.std(0,ddof=1)/np.sqrt(nb)
    a1.plot(tp[1:],e[1:],label=f"{T} K (largo)"); a1.fill_between(tp[1:],(e-es)[1:],(e+es)[1:],alpha=0.25)
    a2.plot(tp[1:],rho[1:],label=f"{T} K")
a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlabel("t (ps)"); a1.set_ylabel("η(t) (mPa·s)"); a1.set_title("η(t), corridas largas (banda: sem)"); a1.legend(fontsize=8)
a2.set_xscale("log"); a2.axhline(0,color="k",lw=0.5); a2.set_xlabel("t (ps)"); a2.set_ylabel("C_P/C_P(0)"); a2.set_title("ACF de presión, corridas largas"); a2.legend(fontsize=8)
# ---- (2) 1000 K, 6 replicas ----
d="cal_1000K_2e10"; T=1000; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); B=[]
for r in range(1,7):
    J=np.loadtxt(f"{d}/cal_r{r}_hflux.dat",comments="#"); dtj=(J[1,0]-J[0,0])*1e-3; m=len(J)//10
    B+=blocks([J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)],int(25/dtj),m)
B=np.array(B); K=(np.cumsum(B,1)-0.5*B[:,:1])*dtj/(V*1e-30*kB*T**2)*1e-12; tj=np.arange(K.shape[1])*dtj
s5=(tj>=0.5)&(tj<=2); i2,i20=int(2/dtj),int(20/dtj)
print("\n===== 1000 K, 6 replicas x 10 bloques de 100 ps")
for grp,idx in (("r1-r3",range(0,30)),("r4-r6",range(30,60)),("r1-r6",range(0,60))):
    Kg=K[list(idx)]; Bg=B[list(idx)]; g=Kg[:,s5].mean(1); Dd=Kg[:,i20]-Kg[:,i2]; tail=Bg[:,i2:i20+1].mean(1)/Bg[:,0].mean()
    print(f"  {grp}: kappa[0.5-2] = {g.mean():.3f} +- {sem(g):.3f} | k(20)-k(2) = {Dd.mean():+.3f} +- {sem(Dd):.3f} | cola C/C0 [2,20] = {1e4*tail.mean():+.2f} +- {1e4*sem(tail):.2f} e-4")
print("  por replica k(20)-k(2): "+" ".join(f"r{r+1}:{(K[r*10:(r+1)*10,i20]-K[r*10:(r+1)*10,i2]).mean():+.2f}" for r in range(6)))
print("  por replica kappa[0.5-2]: "+" ".join(f"r{r+1}:{K[r*10:(r+1)*10][:,s5].mean():.3f}" for r in range(6)))
for nm,idx in (("r1-r3",range(0,30)),("r4-r6",range(30,60))): a3.plot(tj[1:],K[list(idx)].mean(0)[1:],label=nm)
a3.plot(tj[1:],K.mean(0)[1:],"k",lw=2,label="r1-r6"); a3.set_xscale("log"); a3.set_xlim(0.01,25); a3.set_ylim(0,1.6); a3.axvspan(0.5,2,color="gray",alpha=0.12)
a3.set_xlabel("t (ps)"); a3.set_ylabel("κ(t) (W/m/K)"); a3.set_title("κ(t) a 1000 K, 6 réplicas"); a3.legend(fontsize=8)
fig.tight_layout(); fig.savefig("eta_largo.png",dpi=140)
