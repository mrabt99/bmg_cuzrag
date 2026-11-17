"""Fase 3, tanda B: eta y D en el liquido (900-1200 K) para las 3 tasas, con el criterio de eta_T.py.
# COPIA de eta_tasas.py (17/11/2026) con la 1e12 UNIFICADA (cal_*_1e12u) en lugar de la vieja; salidas eta_tasas_1e12u.*. El original no se toca.
Siempre 3 replicas de 1 ns (r1-r3) por caso, bloques de 100 ps, para comparar las tasas en igualdad de condiciones
(2e10 a 1000 K tiene 6 replicas y a 900 K corridas de 3 ns: aca se usan solo r1-r3 de 1 ns).
eta: ACF de pxy,pxz,pyz; t1 = primer tiempo con C(t)/C(0) < 0.02; eta = media de eta(t) en [t1, 3 t1] (si 3 t1 <= TMAX);
     error = sem de los 30 bloques; meseta si |pendiente| < 2*error bootstrap. Ademas eta(t) a 10 y 20 ps fijos.
     Prefactor V/(kB T) con la T media medida en cada replica.
     OJO a 900 K: con 1 ns y bloques de 100 ps eta se subestima (2e10: 37 a 20 ps contra 42.3 con 3 ns x 3, eta_largo.py);
     el sesgo deberia ser el mismo en las tres tasas, asi que sirve para compararlas, no como valor absoluto.
D: pendiente del MSD en la segunda mitad de cada replica, por especie (1e-9 m^2/s). Stokes-Einstein: eta*D_Cu/T.
Uso (desde transporte/): python eta_tasas.py -> eta_tasas.txt, .csv, .png"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
kB=1.380649e-23; NB=10; TMAX=30.0; rng=np.random.default_rng(0); TAGS=["cal_r1","cal_r2","cal_r3"]
RATES=[("2e10",2e10),("1e11",1e11),("1e12u",1e12)]; TS=(900,1000,1100,1200)
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
def caso(T,rate):
    d=f"cal_{T}K_{rate}"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); E=[]; C=[]; D=[]; TM=[]
    for tg in TAGS:
        P=np.loadtxt(f"{d}/{tg}_press.dat",comments="#"); dt=(P[1,0]-P[0,0])*1e-3; nl=int(TMAX/dt)+1
        Tm=np.loadtxt(f"{d}/{tg}_thermo.dat",comments="#")[:,1].mean(); TM.append(Tm)
        cols=[P[:,i]*1e5 for i in (4,5,6)]; m=len(cols[0])//NB
        B=np.array([np.mean([acf(c[b*m:(b+1)*m],nl) for c in cols],axis=0) for b in range(NB)]); C.append(B)
        E.append(V*1e-30/(kB*Tm)*1e-12*(np.cumsum(B,1)-0.5*B[:,:1])*dt*1e3)   # mPa s
        ms=np.loadtxt(f"{d}/{tg}_msd.dat",comments="#"); tm=ms[:,0]*1e-3; h=len(tm)//2
        D.append([np.polyfit(tm[h:],ms[h:,i],1)[0]/6*10 for i in (1,2,3)])
    E=np.concatenate(E); C=np.concatenate(C).mean(0); t=np.arange(E.shape[1])*dt; D=np.array(D)
    below=np.nonzero(C/C[0]<0.02)[0]; t1=t[below[0]] if len(below) else np.nan
    if np.isfinite(t1) and 3*t1<=TMAX:
        s=(t>=t1)&(t<=3*t1); g=E[:,s].mean(1); eta=g.mean(); ee=g.std(ddof=1)/np.sqrt(len(g))
        sl=lambda M: np.polyfit(t[s],M[s],1)[0]; s0=sl(E.mean(0)); se=np.std([sl(E[rng.integers(0,len(E),len(E))].mean(0)) for _ in range(400)])
    else: eta=ee=s0=se=np.nan
    fix=[(E[:,int(round(x/dt))].mean(),E[:,int(round(x/dt))].std(ddof=1)/np.sqrt(len(E))) for x in (10,20)]
    return dict(t=t,E=E.mean(0),Eb=E,C=C/C[0],t1=t1,eta=eta,ee=ee,s0=s0,se=se,fix=fix,D=D.mean(0),De=D.std(0,ddof=1)/np.sqrt(len(D)),rv=dict(RATES)[rate],Tm=np.mean(TM))
R={(T,r):caso(T,r) for T in TS for r,_ in RATES}; out=[]
out.append("eta (mPa s), 3 replicas x 1 ns, bloques de 100 ps")
out.append(" T    tasa  t1(ps)  eta[t1,3t1]  sem    pendiente+-boot   meseta   eta(10ps)     eta(20ps)     | D_Zr   D_Cu   D_Ag (1e-9 m2/s, +-sem rep.)  | eta*D_Cu/T")
for T in TS:
    for r,_ in RATES:
        q=R[T,r]; ok="si" if abs(q["s0"])<2*q["se"] else ("NO" if np.isfinite(q["s0"]) else "--")
        se_=f"{q['eta']:9.2f}  {q['ee']:5.2f}   {q['s0']:+.3f}+-{q['se']:.3f}" if np.isfinite(q["eta"]) else f"{'(sin [t1,3t1] <= '+str(int(TMAX))+' ps)':>38s}"
        out.append(f"{T:5d} {r}  {q['t1']:6.2f}  {se_}   {ok:5s}  "+"  ".join(f"{a:6.2f}+-{b:4.2f}" for a,b in q["fix"])
                   +"  | "+"  ".join(f"{a:.3f}+-{b:.3f}" for a,b in zip(q["D"],q["De"]))+f"  | {(q['eta'] if np.isfinite(q['eta']) else q['fix'][1][0])*q['D'][1]/T:.4f}")
out.append(""); out.append("efecto de la tasa: diferencia relativa contra 2e10 (eta con el criterio [t1,3t1]; a 900 K, eta(20 ps); D_Cu)")
out.append(" T    eta 1e11/2e10      eta 1e12/2e10      D_Cu 1e11/2e10     D_Cu 1e12/2e10")
for T in TS:
    def v(q): return (q["eta"],q["ee"]) if np.isfinite(q["eta"]) else q["fix"][1]
    a=v(R[T,"2e10"]); row=[]
    for r in ("1e11","1e12u"):
        b=v(R[T,r]); x=b[0]/a[0]; ex=x*np.hypot(a[1]/a[0],b[1]/b[0]); row.append(f"{x:.3f}+-{ex:.3f} ({abs(x-1)/ex:.1f}s)")
    da=R[T,"2e10"]["D"][1],R[T,"2e10"]["De"][1]
    for r in ("1e11","1e12u"):
        b=R[T,r]["D"][1],R[T,r]["De"][1]; x=b[0]/da[0]; ex=x*np.hypot(da[1]/da[0],b[1]/b[0]); row.append(f"{x:.3f}+-{ex:.3f} ({abs(x-1)/ex:.1f}s)")
    out.append(f"{T:5d}  "+"   ".join(row))
# correccion a la T nominal: las T medias difieren entre tasas hasta ~4 K (1000 K: 1002.5 en 2e10, 998.7 en 1e11) y D, eta dependen fuerte de T.
# D: ln D_Cu(1/T) cuadratico ajustado a 2e10 (900-1200 K); eta: VFT de eta_vft.py (B = 1132 K, T0 = 586 K).
x=np.array([1/R[T,"2e10"]["Tm"] for T in TS]); p=np.polyfit(x,np.log([R[T,"2e10"]["D"][1] for T in TS]),2)
dlnD=lambda T: -(2*p[0]/T+p[1])/T**2; BV,T0V=1132.0,586.0
out.append(""); out.append("corregido a T nominal: D_Cu*exp(dlnD/dT*(T-Tm)), eta*exp(-B/(T-T0)^2*(T-Tm)); cocientes contra 2e10 (mismos errores relativos)")
out.append(" T    tasa   Tmed     D_Cu corr   eta corr   |  D_Cu/2e10     eta/2e10")
for T in TS:
    def cor(q):
        e=q["eta"] if np.isfinite(q["eta"]) else q["fix"][1][0]; dT=T-q["Tm"]
        return q["D"][1]*np.exp(dlnD(T)*dT), e*np.exp(-BV/(T-T0V)**2*dT)
    a=cor(R[T,"2e10"])
    for r,_ in RATES:
        b=cor(R[T,r]); out.append(f"{T:5d} {r}  {R[T,r]['Tm']:7.1f}   {b[0]:.3f}      {b[1]:6.2f}    |  {b[0]/a[0]:.3f}         {b[1]/a[1]:.3f}")
open("eta_tasas_1e12u.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
def corT(q,T):
    e=q["eta"] if np.isfinite(q["eta"]) else q["fix"][1][0]; dT=T-q["Tm"]
    return q["D"][1]*np.exp(dlnD(T)*dT), e*np.exp(-BV/(T-T0V)**2*dT)
with open("eta_tasas_1e12u.csv","w") as f:
    f.write("T_K,tasa,t1_ps,eta_mPas,sem,pendiente,err_pend,eta10ps,sem10,eta20ps,sem20,D_Zr,D_Cu,D_Ag,semD_Zr,semD_Cu,semD_Ag,T_medida,DCu_Tnom,eta_Tnom\n")
    for (T,r),q in R.items():
        f.write(",".join([str(T),r]+[f"{x:.5g}" for x in (q["t1"],q["eta"],q["ee"],q["s0"],q["se"],*q["fix"][0],*q["fix"][1],*q["D"],*q["De"],q["Tm"],*corT(q,T))])+"\n")
fig,ax=plt.subplots(1,3,figsize=(14,3.8)); cols={"2e10":"C0","1e11":"C1","1e12u":"C2"}; ls={900:"-",1000:"--",1100:"-.",1200:":"}
for (T,r),q in R.items():
    ax[0].plot(q["t"],q["E"],ls[T],color=cols[r],label=f"{r}" if T==900 else None)
    ax[1].plot(q["t"],q["C"],ls[T],color=cols[r])
for r,_ in RATES:
    ax[2].errorbar([1000/T for T in TS],[R[T,r]["eta"] if np.isfinite(R[T,r]["eta"]) else R[T,r]["fix"][1][0] for T in TS],
                   [R[T,r]["ee"] if np.isfinite(R[T,r]["eta"]) else R[T,r]["fix"][1][1] for T in TS],fmt="o-",color=cols[r],capsize=3,label=f"{r} K/s")
ax[0].set_xscale("log"); ax[0].set_xlim(0.05,TMAX); ax[0].set_xlabel("t (ps)"); ax[0].set_ylabel(r"$\eta(t)$ (mPa s)"); ax[0].legend(fontsize=7,title="900 (-), 1000 (--), 1100 (-.), 1200 K (:)",title_fontsize=6)
ax[1].set_xscale("log"); ax[1].set_xlim(0.05,TMAX); ax[1].axhline(0.02,color="k",lw=0.6); ax[1].set_xlabel("t (ps)"); ax[1].set_ylabel("C(t)/C(0) presion")
ax[2].set_yscale("log"); ax[2].set_xlabel("1000/T (1/K)"); ax[2].set_ylabel(r"$\eta$ (mPa s) [900 K: $\eta$(20 ps)]"); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig("eta_tasas_1e12u.png",dpi=150)
