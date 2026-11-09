"""Chequeo de la subida de kappa(t) a tiempos largos (hipotesis 1, 2 y 4). 2e10, T = 300..1200 K.
Para cada T (30 bloques de 100 ps = 3 replicas x 10, ACF hasta 25 ps, trapecio):
 (a) cola de la ACF: media de C(t)/C(0) en [2,20] ps +- sem entre bloques, y fraccion de bloques con media de cola > 0
 (b) subida D = kappa(20) - kappa(2): global +- sem, y por replica
 (c) lo mismo restando la media de J de cada replica (1 ns) antes de correlacionar
 (d) sesgo predicho por la media no restada: pref * <mean_J^2> * 18 ps
Escribe cola_acf.png (ACF cruda 2-20 ps con banda de sem, por T)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
kB=1.380649e-23; eV=1.602176634e-19; NB=10; TAGS=["cal_r1","cal_r2","cal_r3"]; TS=(300,400,500,600,1000,1200)
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
fig,axs=plt.subplots(2,3,figsize=(13,6.5),sharex=True); axs=axs.ravel()
print(" T   | cola C/C0 [2,20] ps (x1e-4)  frac>0 | D=k(20)-k(2) sin restar: global, reps | restando media: global, reps | sesgo media predicho")
for ax,T in zip(axs,TS):
    d=f"cal_{T}K_2e10"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); pref=None
    RES={}
    for sub in (False,True):
        BL=[]; bias=[]
        for tg in TAGS:
            J=np.loadtxt(f"{d}/{tg}_hflux.dat",comments="#"); dt=(J[1,0]-J[0,0])*1e-3; nl=int(25/dt)
            X=[J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)]
            if sub: X=[x-x.mean() for x in X]
            else: bias.append(np.mean([x.mean()**2 for x in X]))
            m=len(X[0])//NB
            BL+=[np.mean([acf(x[b*m:(b+1)*m],nl) for x in X],axis=0) for b in range(NB)]
        BL=np.array(BL); pref=dt/(V*1e-30*kB*T**2)*1e-12
        K=(np.cumsum(BL,1)-0.5*BL[:,:1])*pref; t=np.arange(nl)*dt; i2,i20=int(2/dt),int(20/dt)
        C0=BL[:,0].mean(); tail=BL[:,i2:i20+1].mean(1)/C0
        D=K[:,i20]-K[:,i2]; Dr=[D[r*NB:(r+1)*NB].mean() for r in range(3)]
        RES[sub]=(tail,D,Dr,BL,C0,t)
        if not sub: bpred=np.mean(bias)*18/dt*pref    # pref incluye dt: sesgo = <m^2>*18ps / (V kB T^2)
    tail,D,Dr,BL,C0,t=RES[False]; tailS,DS,DrS,_,_,_=RES[True]
    sem=lambda a: a.std(ddof=1)/np.sqrt(len(a))
    print(f"{T:5d}| {1e4*tail.mean():+6.2f} +- {1e4*sem(tail):5.2f}   {np.mean(tail>0):.2f}  | {D.mean():+.3f} +- {sem(D):.3f}  "+" ".join(f"{x:+.2f}" for x in Dr)
          +f" | {DS.mean():+.3f} +- {sem(DS):.3f}  "+" ".join(f"{x:+.2f}" for x in DrS)+f" | {bpred:+.4f} W/m/K")
    s=(t>=1.5)&(t<=20); mC=BL.mean(0)/C0; sC=BL.std(0,ddof=1)/np.sqrt(len(BL))/C0
    ax.fill_between(t[s],(mC-2*sC)[s],(mC+2*sC)[s],color="C0",alpha=0.25,label="±2 sem (30 bloques)")
    ax.plot(t[s],mC[s],color="C0",lw=0.8,label="C(t)/C(0) media"); ax.axhline(0,color="k",lw=0.6)
    ax.set_title(f"{T} K",fontsize=10); ax.set_ylim(-0.006,0.006)
for ax in axs[3:]: ax.set_xlabel("t (ps)")
axs[0].set_ylabel("C(t)/C(0)"); axs[3].set_ylabel("C(t)/C(0)"); axs[0].legend(fontsize=7)
fig.suptitle("Cola de la ACF del flujo de calor (2e10, sin restar la media)",fontsize=11); fig.tight_layout(); fig.savefig("cola_acf.png",dpi=140)
