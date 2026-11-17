"""Fases 3-4: kappa(T) de 300 a 1200 K para las 3 tasas (tandas A, B y C + corridas previas de 2e10). Criterio:
# COPIA de kappa_tasas.py (17/11/2026) con la 1e12 UNIFICADA (cal_*_1e12u) en lugar de la vieja; salidas kappa_tasas_1e12u.*. El original no se toca.
trapecio, bloques de 100 ps (10 por replica de 1 ns), ventana [0.5,2] ps, test de pendiente (|pendiente| < 2*error bootstrap),
sistematico = std de ventanas [0.5,1],[1,2],[0.5,1.5],[0.75,2].
kappa = V/(3 kB T^2) * integral: se da con T nominal y con la T media medida en cada replica (columna principal: T medida).
Efecto de la tasa: diferencias contra 2e10 y pendiente de kappa contra log10(tasa), por T y combinada por zona
(vidrio 300-600 K, cerca de Tg 700-900 K, liquido 1000-1200 K). Robustez solo para el vidrio. Siempre 3 replicas (r1-r3).
Deriva: T y Ep en los primeros y ultimos 100 ps de la NVE.
Uso (desde transporte/): python kappa_tasas.py -> kappa_tasas.txt, .csv, .png"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
kB=1.380649e-23; eV=1.602176634e-19; kBeV=8.617333e-5; NB=10; NAT=65536; TAGS=["cal_r1","cal_r2","cal_r3"]
W0,W1=0.5,2.0; WINS=[(0.5,1),(1,2),(0.5,1.5),(0.75,2)]
RATES=[("2e10",2e10),("1e11",1e11),("1e12u",1e12)]; TS=(300,400,500,600,700,800,900,1000,1100,1200)
ZONAS=[("vidrio 300-600 K",(300,400,500,600)),("cerca de Tg 700-900 K",(700,800,900)),("liquido 1000-1200 K",(1000,1100,1200))]; GLASS=ZONAS[0][1]
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
def caso(T,rate):
    d=f"cal_{T}K_{rate}"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); K=[]; Tm=[]; dr=[]
    for tg in TAGS:
        J=np.loadtxt(f"{d}/{tg}_hflux.dat",comments="#"); th=np.loadtxt(f"{d}/{tg}_thermo.dat",comments="#")
        dt=(J[1,0]-J[0,0])*1e-3; nl=int(2.5/dt); tm=th[:,1].mean(); Tm.append(tm)
        X=[J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)]; m=len(X[0])//NB
        B=np.array([np.mean([acf(x[b*m:(b+1)*m],nl) for x in X],axis=0) for b in range(NB)])
        K.append((np.cumsum(B,1)-0.5*B[:,:1])*dt/(V*1e-30*kB*T**2)*1e-12)   # con T nominal
        dr.append((th[:100,1].mean(),th[-100:,1].mean(),th[:100,2].mean()/NAT,th[-100:,2].mean()/NAT,th[:,4].max()-th[:,4].min()))
    return K,dt,V,np.array(Tm),np.array(dr)
rng=np.random.default_rng(0); out=[]; R={}
out.append("kappa por trapecio, ventana [0.5,2] ps (W/m/K). Principal: T medida por replica; k_nom: T nominal")
out.append(" T   tasa   <T>     k      sem   sem_rep  sist    k_nom  pend[W/m/K/ps] +-boot  meseta  replicas")
for T in TS:
    for rate,rv in RATES:
        K,dt,V,Tm,dr=caso(T,rate); Kn=np.concatenate(K); Km=np.concatenate([k*(T/tm)**2 for k,tm in zip(K,Tm)])
        t=np.arange(Kn.shape[1])*dt; sel=(t>=W0)&(t<=W1); nb=len(Km)
        g=Km[:,sel].mean(1); k=g.mean(); e=g.std(ddof=1)/np.sqrt(nb); kn=Kn[:,sel].mean()
        rep=g.reshape(3,NB).mean(1); er=rep.std(ddof=1)/np.sqrt(3)
        sist=np.std([Km[:,(t>=a)&(t<=b)].mean() for a,b in WINS],ddof=1)
        sl=lambda M: np.polyfit(t[sel],M[sel],1)[0]
        s0=sl(Km.mean(0)); se=np.std([sl(Km[rng.integers(0,nb,nb)].mean(0)) for _ in range(400)])
        R[T,rate]=dict(k=k,e=e,er=er,sist=sist,kn=kn,s0=s0,se=se,ok=abs(s0)<2*se,Tm=Tm,dr=dr,V=V,rv=rv,t=t,m=Km.mean(0),Km=Km,Kn=Kn,rep=rep)
        out.append(f"{T:4d} {rate} {Tm.mean():7.1f} {k:.3f}  {e:.3f}  {er:.3f}  {sist:.3f}  {kn:.3f}  {s0:+.4f}      {se:.4f}  {'si' if abs(s0)<2*se else 'NO':3s}   "+" ".join(f"{x:.3f}" for x in rep))
out.append(""); out.append("efecto de la tasa (error: sem de bloques; entre parentesis, sigmas)")
out.append(" T    k(1e11)-k(2e10)        k(1e12)-k(2e10)        pendiente dk/dlog10(tasa)")
SL={}; D={}
for T in TS:
    a,b,c=(R[T,r] for r,_ in RATES); d1=b["k"]-a["k"]; e1=np.hypot(a["e"],b["e"]); d2=c["k"]-a["k"]; e2=np.hypot(a["e"],c["e"])
    x=np.log10([r["rv"] for r in (a,b,c)]); y=np.array([r["k"] for r in (a,b,c)]); w=1/np.array([r["e"] for r in (a,b,c)])**2
    xm=np.sum(w*x)/w.sum(); s=np.sum(w*(x-xm)*y)/np.sum(w*(x-xm)**2); es=1/np.sqrt(np.sum(w*(x-xm)**2)); SL[T]=(s,es); D[T]=(d2,e2)
    out.append(f"{T:4d}  {d1:+.3f} +- {e1:.3f} ({abs(d1)/e1:.1f})    {d2:+.3f} +- {e2:.3f} ({abs(d2)/e2:.1f})    {s:+.4f} +- {es:.4f} ({abs(s)/es:.1f})")
for zn,zt in ZONAS:
    s,es=map(np.array,zip(*[SL[T] for T in zt])); w=1/es**2; S=np.sum(w*s)/w.sum(); ES=1/np.sqrt(w.sum()); chi=np.sum(w*(s-S)**2)
    d,ed=map(np.array,zip(*[D[T] for T in zt])); w2=1/ed**2; DM=np.sum(w2*d)/w2.sum(); EDM=1/np.sqrt(w2.sum())
    out.append(f"combinado {zn:22s}: pendiente = {S:+.4f} +- {ES:.4f} W/m/K por decada ({abs(S)/ES:.1f} sigma), chi2 entre T = {chi:.1f}/{len(s)-1};"
               f"  k(1e12)-k(2e10) = {DM:+.3f} +- {EDM:.3f} ({abs(DM)/EDM:.1f} sigma)")
POOL=np.sqrt(np.mean([np.var(R[T,r]["rep"],ddof=1) for T in GLASS for r,_ in RATES])/3)   # sem de una media de 3 replicas, varianza agrupada (24 gdl)
SEMB=np.mean([R[T,r]["e"] for T in GLASS for r,_ in RATES])
def comb(temps=GLASS,win=(W0,W1),key="Km",err=None):
    s=[];es=[];d=[];ed=[]
    for T in temps:
        v=[]
        for rate,_ in RATES:
            r=R[T,rate]; sel=(r["t"]>=win[0])&(r["t"]<=win[1]); g=r[key][:,sel].mean(1); v.append((g.mean(),err if err else g.std(ddof=1)/np.sqrt(len(g))))
        y,e=map(np.array,zip(*v)); x=np.log10([rv for _,rv in RATES]); w=1/e**2; xm=np.sum(w*x)/w.sum()
        s.append(np.sum(w*(x-xm)*y)/np.sum(w*(x-xm)**2)); es.append(1/np.sqrt(np.sum(w*(x-xm)**2))); d.append(y[2]-y[0]); ed.append(np.hypot(e[0],e[2]))
    s,es,d,ed=map(np.array,(s,es,d,ed)); w=1/es**2; w2=1/ed**2
    return np.sum(w*s)/w.sum(),1/np.sqrt(w.sum()),np.sum(w2*d)/w2.sum(),1/np.sqrt(w2.sum())
out.append(""); out.append(f"robustez del efecto combinado en el vidrio (sem medio de bloques {SEMB:.4f}; sem de 3 replicas con varianza agrupada {POOL:.4f})")
out.append("variante                         pendiente/decada            k(1e12)-k(2e10)")
for lab,kw in (("base (T medida, 0.5-2 ps)",{}),("T nominal",dict(key="Kn")),("sin 600 K",dict(temps=(300,400,500))),("ventana 0.5-1 ps",dict(win=(0.5,1))),("ventana 1-2 ps",dict(win=(1,2))),("error entre replicas",dict(err=POOL)),("sin 600 K, error entre replicas",dict(temps=(300,400,500),err=POOL))):
    a1,a2,a3,a4=comb(**kw); out.append(f"{lab:32s} {a1:+.4f} +- {a2:.4f} ({abs(a1)/a2:.1f})   {a3:+.3f} +- {a4:.3f} ({abs(a3)/a4:.1f})")
out.append(""); out.append("deriva en NVE (media de las 3 replicas): primeros -> ultimos 100 ps. Conservacion de energia: dT = -dEp/(1.5 kB). Energia inherente: dE_IS = 2 dEp")
out.append(" T   tasa   T ini -> fin      dT     dEp (meV/at)   dT = -dEp/(1.5kB)   Ep/at (eV)   V/at (A^3)   rango Etot (eV)")
for T in TS:
    for rate,_ in RATES:
        r=R[T,rate]; q=r["dr"].mean(0); dE=(q[3]-q[2])*1e3
        out.append(f"{T:4d} {rate}  {q[0]:7.2f} -> {q[1]:7.2f}  {q[1]-q[0]:+.2f}   {dE:+.3f}        {-dE*1e-3/(1.5*kBeV):+.2f}             {q[2]:.4f}    {r['V']/NAT:.3f}     {q[4]:.1f}")
open("kappa_tasas_1e12u.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
with open("kappa_tasas_1e12u.csv","w") as f:
    f.write("T_K,tasa,T_medida,kappa,sem,sem_replicas,sistematico,kappa_T_nominal,pendiente,err_pendiente,meseta\n")
    for (T,rate),r in R.items(): f.write(f"{T},{rate},{r['Tm'].mean():.2f},{r['k']:.4f},{r['e']:.4f},{r['er']:.4f},{r['sist']:.4f},{r['kn']:.4f},{r['s0']:+.5f},{r['se']:.5f},{'si' if r['ok'] else 'no'}\n")
fig,ax=plt.subplots(1,2,figsize=(11,3.8)); off={"2e10":-8,"1e11":0,"1e12u":8}
for rate,_ in RATES:
    ax[0].errorbar([R[T,rate]["Tm"].mean()+off[rate] for T in TS],[R[T,rate]["k"] for T in TS],[np.hypot(R[T,rate]["e"],R[T,rate]["sist"]) for T in TS],fmt="o-",capsize=3,label=f"{rate} K/s")
    r=R[500,rate]; ax[1].plot(r["t"],r["m"],label=f"{rate} K/s")
ax[0].axvspan(664,692,color="0.88",label="Tg (664-692 K)"); ax[0].set_xlabel("T (K)"); ax[0].set_ylabel(r"$\kappa$ (W/m/K)"); ax[0].set_ylim(0.6,1.0); ax[0].legend(fontsize=8)
ax[1].axvspan(W0,W1,color="0.9"); ax[1].set_xscale("log"); ax[1].set_xlim(0.02,2.5); ax[1].set_xlabel("t (ps)"); ax[1].set_ylabel(r"$\kappa(t)$ a 500 K (W/m/K)"); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig("kappa_tasas_1e12u.png",dpi=150)
# --- Fase 4: cambio de pendiente cerca de Tg (todas las tasas juntas, T medida, error = sem de bloques) ---
Tq=np.array([R[T,r]["Tm"].mean() for T in TS for r,_ in RATES]); kq=np.array([R[T,r]["k"] for T in TS for r,_ in RATES]); wq=1/np.array([R[T,r]["e"] for T in TS for r,_ in RATES])**2
o2=["","cambio de pendiente (30 puntos: 10 T x 3 tasas)"]
for T in TS:
    w=np.array([1/R[T,r]["e"]**2 for r,_ in RATES]); v=np.array([R[T,r]["k"] for r,_ in RATES]); o2.append(f"  {T:5d} K: kappa medio de las 3 tasas = {np.sum(w*v)/w.sum():.3f} +- {1/np.sqrt(w.sum()):.3f}")
A=np.vstack([np.ones_like(Tq),Tq]).T; c=np.linalg.lstsq(A*np.sqrt(wq)[:,None],kq*np.sqrt(wq),rcond=None)[0]; chi1=np.sum(wq*(kq-A@c)**2)
o2.append(f"recta unica: chi2 = {chi1:.1f}/{len(kq)-2}")
BR=[]
for Tb in np.arange(500,1151,10):
    X=np.vstack([np.ones_like(Tq),Tq-Tb,np.clip(Tq-Tb,0,None)]).T; cb=np.linalg.lstsq(X*np.sqrt(wq)[:,None],kq*np.sqrt(wq),rcond=None)[0]; BR.append((np.sum(wq*(kq-X@cb)**2),Tb,cb))
chi0,Tb0,cb0=min(BR,key=lambda z:z[0]); okb=[b[1] for b in BR if b[0]<=chi0+1]
o2.append(f"recta quebrada continua: Tb = {Tb0} K (chi2+1: {min(okb)}-{max(okb)} K), chi2 = {chi0:.1f}/{len(kq)-4}; pendiente {100*cb0[1]:+.4f} antes y {100*(cb0[1]+cb0[2]):+.4f} W/m/K por 100 K despues")
o2.append("chi2 con quiebre forzado: "+", ".join(f"{b[1]} K {b[0]:.1f}" for b in BR if b[1] in (680,800,900,1000)))
open("kappa_tasas_1e12u.txt","a").write("\n".join(o2)+"\n"); print("\n".join(o2))
