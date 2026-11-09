"""eta(T), kappa(T) y D(T) en el liquido a 2e10 con criterio de ventana objetivo (corridas de 1 ns, bloques de 100 ps).
eta: ACF de pxy,pxz,pyz promediada (30+ bloques). t1 = primer tiempo con C(t)/C(0) < 0.02 (ACF media).
     eta = promedio de eta(t) en [t1, 3*t1]; error = sem entre bloques; meseta si |pendiente| < 2*error bootstrap.
     Control: regla tipo Zhang (corte donde std/media entre bloques = 40 %, valor = media en [tc/2, tc]).
kappa: trapecio, ventana [0.5, 2] ps (como mail_seba/kappa_vs_T.py). D: pendiente del MSD en la 2da mitad.
Uso (desde transporte/): python eta_T.py -> eta_T.txt, eta_T.csv, eta_T.png. 900 K sale de eta_largo.py (500 ps)."""
import glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
kB=1.380649e-23; eV=1.602176634e-19; NB=10; TMAX=20.0; rng=np.random.default_rng(0)
TS=[1000,1100,1200,1300,1400,1500,1600]
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
def blocks(cols,nl):
    m=len(cols[0])//NB; return np.array([np.mean([acf(c[b*m:(b+1)*m],nl) for c in cols],axis=0) for b in range(NB)])
def ventana(E,t,a,b):
    s=(t>=a)&(t<=b); g=E[:,s].mean(1); v=g.mean(); e=g.std(ddof=1)/np.sqrt(len(g))
    sl=lambda M: np.polyfit(t[s],M[s],1)[0]; s0=sl(E.mean(0))
    se=np.std([sl(E[rng.integers(0,len(E),len(E))].mean(0)) for _ in range(400)]); return v,e,s0,se
rows=[]; out=[]; fig,ax=plt.subplots(1,3,figsize=(14,3.8))
out.append(" T    nrep  t1(ps)  eta[t1,3t1]   sem   pend+-boot     meseta  Zhang40  |  kappa   sem  meseta |  D_Zr  D_Cu  D_Ag (1e-9 m2/s)")
for T in TS:
    d=f"cal_{T}K_2e10"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); tags=sorted(glob.glob(f"{d}/cal_r*_press.dat"))
    E=[];K=[];C=[];D=[]
    for f in tags:
        tg=f[:-len("_press.dat")]
        P=np.loadtxt(f,comments="#"); dtp=(P[1,0]-P[0,0])*1e-3; nl=int(TMAX/dtp)+1
        B=blocks([P[:,i]*1e5 for i in (4,5,6)],nl); C.append(B)
        E.append(V*1e-30/(kB*T)*1e-12*(np.cumsum(B,1)-0.5*B[:,:1])*dtp*1e3)
        J=np.loadtxt(f"{tg}_hflux.dat",comments="#"); dtj=(J[1,0]-J[0,0])*1e-3
        BJ=blocks([J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)],int(2.5/dtj))
        K.append((np.cumsum(BJ,1)-0.5*BJ[:,:1])*dtj/(V*1e-30*kB*T**2)*1e-12)
        m=np.loadtxt(f"{tg}_msd.dat",comments="#"); tm=m[:,0]*1e-3; h=len(tm)//2
        D.append([np.polyfit(tm[h:],m[h:,i],1)[0]/6*10 for i in (1,2,3)])
    E=np.concatenate(E); K=np.concatenate(K); C=np.concatenate(C).mean(0); D=np.mean(D,0)
    t=np.arange(E.shape[1])*dtp; tk=np.arange(K.shape[1])*dtj
    t1=t[np.argmax(C/C[0]<0.02)]
    eta,ee,s0,se=ventana(E,t,t1,3*t1)
    rel=E.std(0,ddof=1)/np.abs(E.mean(0)); i40=np.argmax((rel>0.40)&(t>0.2)); tc=t[i40] if i40>0 else t[-1]
    zh=E[:,(t>=tc/2)&(t<=tc)].mean()
    k,ke,ks0,kse=ventana(K,tk,0.5,2.0)
    rows.append((T,len(tags),t1,eta,ee,s0,se,zh,k,ke,*D))
    out.append(f"{T:5d} {len(tags):3d}  {t1:6.2f}  {eta:9.2f}  {ee:6.2f}  {s0:+7.3f}+-{se:.3f}  {'si' if abs(s0)<2*se else 'NO':5s} {zh:7.2f}  |  {k:.3f} {ke:.3f}  {'si' if abs(ks0)<2*kse else 'NO':4s} | "
               +"  ".join(f"{x:.3f}" for x in D))
    l,=ax[0].plot(t,E.mean(0),label=f"{T} K"); ax[0].axvspan(t1,3*t1,color=l.get_color(),alpha=0.12)
    ax[1].plot(t,C/C[0],color=l.get_color())
out.append("900 K (eta_largo.py, 3x3 ns, bloques de 500 ps, ventana [20,60] ps): eta = 42.3 +- 2.7 mPa s, kappa = 0.824 +- 0.009")
R=np.array([r for r in rows],float)
ax[2].errorbar(1000/R[:,0],R[:,3],R[:,4],fmt="o",label="1 ns x 3-6 rep.")
ax[2].errorbar([1000/900],[42.3],[2.7],fmt="s",label="900 K (3 ns x 3)")
x=1000/R[:,0]; p=np.polyfit(x,np.log(R[:,3]),1); Ea=p[0]*1000*8.617333e-5
ax[2].plot(x,np.exp(np.polyval(p,x)),"--",lw=0.8,label=f"Arrhenius 1000-1600 K: Ea={Ea:.2f} eV")
out.append(f"Arrhenius (1000-1600 K): Ea = {Ea:.3f} eV; extrapolado a 900 K: {np.exp(np.polyval(p,1000/900)):.1f} mPa s (medido 42.3)")
for Ti,DCu,et in zip(R[:,0],R[:,11],R[:,3]):
    out.append(f"Stokes-Einstein T={Ti:.0f}: eta*D_Cu/T = {et*DCu/Ti:.4f} (mPa s 1e-9 m2/s / K)")
ax[0].set_xscale("log"); ax[0].set_xlabel("t (ps)"); ax[0].set_ylabel(r"$\eta(t)$ (mPa s)"); ax[0].legend(fontsize=7)
ax[1].set_xscale("log"); ax[1].set_xlabel("t (ps)"); ax[1].set_ylabel("C(t)/C(0) presion"); ax[1].axhline(0.02,color="k",lw=0.6)
ax[2].set_yscale("log"); ax[2].set_xlabel("1000/T (1/K)"); ax[2].set_ylabel(r"$\eta$ (mPa s)"); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig("eta_T.png",dpi=150)
with open("eta_T.csv","w") as f:
    f.write("T_K,nrep,t1_ps,eta_mPas,sem,pendiente,err_pend,eta_zhang40,kappa,kappa_sem,D_Zr,D_Cu,D_Ag\n")
    for r in rows: f.write(",".join(f"{v:.5g}" for v in r)+"\n")
open("eta_T.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
