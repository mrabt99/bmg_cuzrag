"""Curvas GK para 2e10 en todas las T (300-1200 K), con enfasis en 700/800/900 K. Trapecio, 30 bloques de 100 ps.
Figura curvas_T.png: (a) kappa(t)  (b) ACF de J normalizada  (c) eta(t)  (d) ACF de presion normalizada
Imprime: chequeos de equilibrado (NPT, NVE), kappa (ventana 0.5-2 ps + test de meseta), eta(t), D por especie."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
kB=1.380649e-23; eV=1.602176634e-19; NB=10; TAGS=["cal_r1","cal_r2","cal_r3"]; N=65536
TS=(300,400,500,600,700,800,900,1000,1200); NEW=(700,800,900); rng=np.random.default_rng(0)
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
def blocks(cols,nl):
    m=len(cols[0])//NB; return [np.mean([acf(c[b*m:(b+1)*m],nl) for c in cols],axis=0) for b in range(NB)]
cm=plt.cm.coolwarm(np.linspace(0,1,len(TS)))
fig,ax=plt.subplots(2,2,figsize=(12,8.5)); a,b,c,d=ax.ravel()
print(" T  | NPT 2a mitad dV/dt(%/ps) dEp/dt(meV/at/ps) | NVE: Tmedia, dEp(fin-ini) meV/at | kappa[0.5-2] sem  pend(+-boot) meseta | eta(t=1,3,5,10,20 ps) mPa s | D(Zr,Cu,Ag) 1e-9 m2/s")
for col,T in zip(cm,TS):
    dd=f"cal_{T}K_2e10"; V=float(open(f"{dd}/cal_vmean.inc").read().split()[3])
    L=np.loadtxt(f"{dd}/cal_npt.log.dat",comments="#"); h=len(L)//2; tt=(L[:,0]-L[0,0])*1e-3
    sv=np.polyfit(tt[h:],L[h:,3],1)[0]/L[h:,3].mean()*100; se=np.polyfit(tt[h:],L[h:,4]/N,1)[0]*1e3
    BJ=[];BP=[];th=[];Ds=[]
    for tg in TAGS:
        J=np.loadtxt(f"{dd}/{tg}_hflux.dat",comments="#"); dtj=(J[1,0]-J[0,0])*1e-3
        BJ+=blocks([J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)],int(3/dtj))
        P=np.loadtxt(f"{dd}/{tg}_press.dat",comments="#"); dtp=(P[1,0]-P[0,0])*1e-3
        BP+=blocks([P[:,i]*1e5 for i in (4,5,6)],int(25/dtp))
        x=np.loadtxt(f"{dd}/{tg}_thermo.dat",comments="#"); k=len(x)//10; th.append((x[:,1].mean(),(x[-k:,2].mean()-x[:k,2].mean())/N*1e3))
        m=np.loadtxt(f"{dd}/{tg}_msd.dat",comments="#"); t=m[:,0]*1e-3; hh=len(t)//2
        Ds.append([np.polyfit(t[hh:],m[hh:,i],1)[0]/6*10 for i in (1,2,3)])
    BJ=np.array(BJ); BP=np.array(BP)
    K=(np.cumsum(BJ,1)-0.5*BJ[:,:1])*dtj/(V*1e-30*kB*T**2)*1e-12; tj=np.arange(K.shape[1])*dtj
    E=(np.cumsum(BP,1)-0.5*BP[:,:1])*dtp*V*1e-30/(kB*T)*1e-12*1e3; tp=np.arange(E.shape[1])*dtp
    sel=(tj>=0.5)&(tj<=2); g=K[:,sel].mean(1); kap=g.mean(); ke=g.std(ddof=1)/np.sqrt(len(g))
    sl=lambda M: np.polyfit(tj[sel],M[sel],1)[0]; s0=sl(K.mean(0)); sb=np.std([sl(K[rng.integers(0,30,30)].mean(0)) for _ in range(300)])
    eta=" ".join(f"{E[:,int(x/dtp)].mean():.1f}" for x in (1,3,5,10,20))
    D=np.mean(Ds,axis=0)
    mark="*" if T in NEW else " "
    print(f"{T:4d}{mark}| {sv:+.5f} {se:+.4f} | {np.mean([x[0] for x in th]):.1f}, "+",".join(f"{x[1]:+.2f}" for x in th)
          +f" | {kap:.3f} {ke:.3f} {s0:+.3f}({sb:.3f}) {'si' if abs(s0)<2*sb else 'NO'} | {eta} | "+" ".join(f"{x:.3f}" for x in D))
    lw=2.0 if T in NEW else 1.0
    a.plot(tj[1:],K.mean(0)[1:],color=col,lw=lw,label=f"{T} K")
    b.plot(tj,BJ.mean(0)/BJ[:,0].mean(),color=col,lw=lw,label=f"{T} K")
    if T>=700:
        e=E.mean(0); es=E.std(0,ddof=1)/np.sqrt(30)
        c.plot(tp[1:],e[1:],color=col,lw=lw,label=f"{T} K"); c.fill_between(tp[1:],(e-es)[1:],(e+es)[1:],color=col,alpha=0.2)
    d.plot(tp[1:],(BP.mean(0)/BP[:,0].mean())[1:],color=col,lw=lw,label=f"{T} K")
a.set_xscale("log"); a.set_xlim(0.005,3); a.set_ylim(0,2.7); a.set_xlabel("t (ps)"); a.set_ylabel("κ(t) (W/m/K)"); a.set_title("κ(t) acumulada (trapecio)"); a.axvspan(0.5,2,color="gray",alpha=0.12); a.legend(fontsize=7,ncol=2)
b.set_xscale("log"); b.set_xlim(0.005,3); b.axhline(0,color="k",lw=0.5); b.set_xlabel("t (ps)"); b.set_ylabel("C_J(t)/C_J(0)"); b.set_title("ACF del flujo de calor")
c.set_xscale("log"); c.set_yscale("log"); c.set_xlim(0.01,25); c.set_xlabel("t (ps)"); c.set_ylabel("η(t) (mPa·s)"); c.set_title("η(t) acumulada, 700-1200 K (banda: sem 30 bloques)"); c.legend(fontsize=7)
d.set_xscale("log"); d.set_xlim(0.002,25); d.axhline(0,color="k",lw=0.5); d.set_xlabel("t (ps)"); d.set_ylabel("C_P(t)/C_P(0)"); d.set_title("ACF de presión (pxy, pxz, pyz)")
fig.suptitle("Green-Kubo, 2e10 K/s, 65536 átomos (700/800/900 K en trazo grueso)",fontsize=11); fig.tight_layout(); fig.savefig("curvas_T.png",dpi=140)
