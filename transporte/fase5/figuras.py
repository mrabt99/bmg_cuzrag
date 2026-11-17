"""Fase 5: figuras finales de transporte (kappa, eta, D) para la tesis.
Lee los resultados ya calculados (no recalcula nada salvo el pico de kappa(t), que necesita las series de J):
  ../kappa_tasas.csv   kappa(T) 300-1200 K, 3 tasas, T medida (kappa_tasas.py)
  ../eta_T.csv         2e10 en 1000-1600 K: eta, kappa (T nominal), D (eta_T.py)
  ../eta_tasas.csv     eta y D_Cu 900-1200 K, 3 tasas, corregidos a T nominal (eta_tasas.py)
  ../eta_largo.log     900 K (3 x 3 ns) y 800 K (3 x 5 ns), 2e10 (eta_largo.py)
  ../eta_vft.txt       parametros VFT (eta_vft.py)
  ../../structural_analysis/Tg_resultados.csv   Tg por tasa (Etapa 1, media de Tg(V) y Tg(Ep))
Escritos a mano (con su fuente): kappa NEMD (nemd/README.md, 2026-11-10) y Mendez et al., JNCS 693 (2026) 124394, Tabla 3.
Salidas: fig_kappa_T.{pdf,png}, fig_eta_D_T.{pdf,png}
Uso (desde transporte/fase5/): python figuras.py"""
import csv, re, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.labelsize":9,"legend.fontsize":7,"xtick.direction":"in","ytick.direction":"in","xtick.top":True,"ytick.right":True})
kB=1.380649e-23; eV=1.602176634e-19; NB=10
RATES=["2e10","1e11","1e12"]; LAB={"2e10":r"$2\times10^{10}$ K/s","1e11":r"$10^{11}$ K/s","1e12":r"$10^{12}$ K/s"}
COL={"2e10":plt.cm.viridis(0.05),"1e11":plt.cm.viridis(0.45),"1e12":plt.cm.viridis(0.75)}; MK={"2e10":"o","1e11":"s","1e12":"^"}; OFF={"2e10":-7,"1e11":0,"1e12":7}
rd=lambda f: list(csv.DictReader(open(f)))
KT=rd("../kappa_tasas.csv"); ET=rd("../eta_T.csv"); EA=rd("../eta_tasas.csv")
TG={r["tasa"]:(float(r["Tg_V_K"])+float(r["Tg_Ep_K"]))/2 for r in rd("../../structural_analysis/Tg_resultados.csv")}
vft=open("../eta_vft.txt").readline(); ETA0,BV,T0V=(float(re.search(p,vft).group(1)) for p in (r"eta0 = ([\d.]+)",r"B = D\*T0 = ([\d.]+)",r", T0 = ([\d.]+) K"))
LOG=open("../eta_largo.log").read()
def largo(T):   # bloque de eta_largo.log para la T dada
    b=LOG.split(f"===== {T} K")[1].split("=====")[0]
    D=[float(x) for x in re.search(r"D\(Zr,Cu,Ag\) = ([\d.]+) ([\d.]+) ([\d.]+)",b).groups()]
    return b,D
b900,D900=largo(900); E900=[float(x) for x in re.search(r"ventana \[20,60\] ps: eta = ([\d.]+) \+- ([\d.]+)",b900).groups()]
K900=[float(x) for x in re.search(r"kappa\[0.5-2 ps\] = ([\d.]+) \+- ([\d.]+)",b900).groups()]
b800,D800=largo(800); E800=[float(x) for x in re.search(r"t=199:([\d.]+)\+-([\d.]+)",b800).groups()]   # eta(200 ps): cota inferior
NEMD=(300,0.730,0.016)   # nemd/README.md: kappa_inf, vidrio de 8192 at. a 1e12, 300 K
MENDEZ=[(300,2.43,0.71),(400,1.37,0.40),(500,1.12,0.59),(600,0.77,0.21),(1000,0.94,np.nan)]   # Tabla 3 y texto (1000 K sin error)

# --- pico de kappa(t) en 2e10 (mismo calculo que kappa_tasas.py, solo hasta 0.3 ps) ---
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
def pico(T):
    d=f"../cal_{T}K_2e10"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); K=[]
    for r in (1,2,3):
        J=np.loadtxt(f"{d}/cal_r{r}_hflux.dat",comments="#"); Tm=np.loadtxt(f"{d}/cal_r{r}_thermo.dat",comments="#")[:,1].mean()
        dt=(J[1,0]-J[0,0])*1e-3; nl=int(0.3/dt); X=[J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)]; m=len(X[0])//NB
        B=np.array([np.mean([acf(x[b*m:(b+1)*m],nl) for x in X],0) for b in range(NB)])
        K.append((np.cumsum(B,1)-0.5*B[:,:1])*dt/(V*1e-30*kB*Tm**2)*1e-12)
    K=np.concatenate(K); i=np.argmax(K.mean(0)); return K[:,i].mean(),K[:,i].std(ddof=1)/np.sqrt(len(K)),i*dt
PICO={T:pico(T) for T in (300,400,500,600,1000)}

# --- recta quebrada (igual que al final de kappa_tasas.py) ---
Tq=np.array([float(r["T_medida"]) for r in KT]); kq=np.array([float(r["kappa"]) for r in KT]); wq=1/np.array([float(r["sem"]) for r in KT])**2
def quebrada(Tb):
    X=np.vstack([np.ones_like(Tq),Tq-Tb,np.clip(Tq-Tb,0,None)]).T; c=np.linalg.lstsq(X*np.sqrt(wq)[:,None],kq*np.sqrt(wq),rcond=None)[0]; return np.sum(wq*(kq-X@c)**2),c
BR=[(quebrada(Tb),Tb) for Tb in np.arange(560,900,2)]; (chi0,cb),Tb0=min(BR,key=lambda z:z[0][0])
IC=[Tb for (c2,_),Tb in BR if c2<=chi0+1]

# ================= Figura 1: kappa(T) =================
fig,ax=plt.subplots(1,2,figsize=(7.2,3.0),gridspec_kw=dict(width_ratios=[1.35,1]))
a=ax[0]
for r in RATES:
    q=[x for x in KT if x["tasa"]==r]; T=np.array([float(x["T_medida"]) for x in q])+OFF[r]
    k=np.array([float(x["kappa"]) for x in q]); e=np.hypot([float(x["sem"]) for x in q],[float(x["sistematico"]) for x in q])
    if r=="2e10":   # 1300-1600 K solo en 2e10 (eta_T.csv: T nominal, error = sem)
        h=[x for x in ET if int(x["T_K"])>=1300]; T=np.r_[T,[float(x["T_K"])+OFF[r] for x in h]]
        k=np.r_[k,[float(x["kappa"]) for x in h]]; e=np.r_[e,[float(x["kappa_sem"]) for x in h]]
    a.errorbar(T,k,e,fmt=MK[r]+"-",ms=3.5,lw=0.8,capsize=2,color=COL[r],label=LAB[r])
    a.axvline(TG[r],color=COL[r],lw=0.7,ls=":")
xx=np.linspace(300,1200,200); a.plot(xx,cb[0]+cb[1]*(xx-Tb0)+cb[2]*np.clip(xx-Tb0,0,None),"k--",lw=0.9,label=f"recta quebrada ($T_b$ = {Tb0:.0f} K)")
a.axvspan(min(IC),max(IC),color="0.9",zorder=0)
a.errorbar([NEMD[0]-14],[NEMD[1]],[NEMD[2]],fmt="*",ms=7,color="C3",capsize=2,label="NEMD ($L\\to\\infty$)")
a.text(TG["2e10"]-8,0.705,r"$T_g$",fontsize=7,ha="right")
a.set_xlabel("$T$ (K)"); a.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)"); a.set_ylim(0.68,0.96); a.set_xlim(250,1650); a.legend(loc="upper left",ncol=1,frameon=False)
a.text(0.97,0.04,"(a)",transform=a.transAxes,ha="right")
a=ax[1]
q=[x for x in KT if x["tasa"]=="2e10" and int(x["T_K"]) in PICO]
a.errorbar([float(x["T_K"]) for x in q],[float(x["kappa"]) for x in q],[np.hypot(float(x["sem"]),float(x["sistematico"])) for x in q],fmt="o-",ms=3.5,lw=0.8,capsize=2,color=COL["2e10"],label="GK, meseta 0.5–2 ps")
a.errorbar(list(PICO),[v[0] for v in PICO.values()],[v[1] for v in PICO.values()],fmt="v--",ms=4,lw=0.8,capsize=2,color="0.35",label=r"GK, máximo de $\kappa(t)$")
a.errorbar([m[0]+10 for m in MENDEZ],[m[1] for m in MENDEZ],[m[2] for m in MENDEZ],fmt="s",mfc="none",ms=4,capsize=2,color="C3",label="Méndez et al. (2026)")
a.errorbar([NEMD[0]-14],[NEMD[1]],[NEMD[2]],fmt="*",ms=7,color="C3",capsize=2,label="NEMD")
a.set_xlabel("$T$ (K)"); a.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)"); a.set_ylim(0,3.3); a.legend(frameon=False,loc="upper right")
a.text(0.04,0.04,"(b)",transform=a.transAxes)
fig.tight_layout(); [fig.savefig(f"fig_kappa_T.{x}",dpi=300) for x in ("pdf","png")]

# ================= Figura 2: eta, D, Stokes-Einstein =================
E2=[(float(x["T_K"]),float(x["eta_mPas"]),float(x["sem"]),float(x["D_Zr"]),float(x["D_Cu"]),float(x["D_Ag"])) for x in ET]
E2=sorted(E2+[(900.0,E900[0],E900[1],*D900)]); E2=np.array(E2)
fig,ax=plt.subplots(1,3,figsize=(7.2,2.7))
a=ax[0]; x=1000/E2[:,0]
a.errorbar(x,E2[:,1],E2[:,2],fmt="o",ms=3.5,capsize=2,color=COL["2e10"],label=LAB["2e10"])
for r in ("1e11","1e12"):
    q=[y for y in EA if y["tasa"]==r and int(y["T_K"])>=1000]
    a.errorbar([1000/float(y["T_K"]) for y in q],[float(y["eta_Tnom"]) for y in q],[float(y["sem"]) for y in q],fmt=MK[r],mfc="none",ms=4,capsize=2,color=COL[r],label=LAB[r])
for r in RATES:   # 900 K con 1 ns: eta(20 ps), subestimado
    y=[y for y in EA if y["tasa"]==r and y["T_K"]=="900"][0]; a.plot(1000/900+0.004*(RATES.index(r)-1),float(y["eta20ps"]),MK[r],mfc="none",ms=3,color=COL[r],alpha=0.6)
a.errorbar([1000/800],[E800[0]],yerr=[[0],[E800[0]*0.6]],lolims=[True],fmt="o",ms=3.5,color=COL["2e10"],capsize=0)
tt=np.linspace(780,1650,200); a.plot(1000/tt,ETA0*np.exp(BV/(tt-T0V)),"k-",lw=0.8,label=f"VFT ($T_0$ = {T0V:.0f} K)")
m=E2[:,0]>=1000; p=np.polyfit(1000/E2[m,0],np.log(E2[m,1]),1); tt=np.linspace(780,1650,50)
a.plot(1000/tt,np.exp(np.polyval(p,1000/tt)),"k:",lw=0.8,label=f"Arrhenius ({p[0]*1000*8.617333e-5:.2f} eV)")
a.set_yscale("log"); a.set_ylim(2,800); a.set_xlabel(r"$1000/T$ (K$^{-1}$)"); a.set_ylabel(r"$\eta$ (mPa s)"); a.legend(frameon=False,loc="upper left",fontsize=6)
a.text(0.96,0.04,"(a)",transform=a.transAxes,ha="right")
a=ax[1]; E3=np.vstack([E2,[800.0,np.nan,np.nan,*D800]]); E3=E3[np.argsort(E3[:,0])]; x=1000/E3[:,0]
for j,(sp,mk) in enumerate((("Zr","o"),("Cu","s"),("Ag","^"))):
    a.plot(x,E3[:,3+j],mk+"-",ms=3.5,lw=0.8,color=f"C{j}",label=f"{sp} (2e10)")
for r in ("1e11","1e12"):
    q=[y for y in EA if y["tasa"]==r]; a.plot([1000/float(y["T_K"]) for y in q],[float(y["DCu_Tnom"]) for y in q],MK[r],mfc="none",ms=4.5,color=COL[r],label=f"Cu ({r})")
a.set_yscale("log"); a.set_ylim(0.02,8); a.set_xlabel(r"$1000/T$ (K$^{-1}$)"); a.set_ylabel(r"$D$ ($10^{-9}$ m$^2$ s$^{-1}$)"); a.legend(frameon=False,fontsize=6,loc="lower left")
a.text(0.96,0.92,"(b)",transform=a.transAxes,ha="right")
a=ax[2]; se=E2[:,1]*E2[:,4]/E2[:,0]; ref=se[E2[:,0]>=1000].mean()
a.errorbar(E2[:,0],se/ref,se/ref*E2[:,2]/E2[:,1],fmt="o-",ms=3.5,lw=0.8,capsize=2,color=COL["2e10"],label=LAB["2e10"])
for r in ("1e11","1e12"):
    q=[y for y in EA if y["tasa"]==r and int(y["T_K"])>=1000]
    v=np.array([float(y["eta_Tnom"])*float(y["DCu_Tnom"])/float(y["T_K"]) for y in q])
    a.errorbar([float(y["T_K"])+OFF[r] for y in q],v/ref,v/ref*np.array([float(y["sem"])/float(y["eta_Tnom"]) for y in q]),fmt=MK[r],mfc="none",ms=4,capsize=2,color=COL[r],label=LAB[r])
a.axhline(1,color="k",lw=0.6); a.set_xlabel("$T$ (K)"); a.set_ylabel(r"$\eta D_{\rm Cu}/T$ (relativo)"); a.legend(frameon=False,fontsize=6)
a.text(0.96,0.92,"(c)",transform=a.transAxes,ha="right")
fig.tight_layout(); [fig.savefig(f"fig_eta_D_T.{x}",dpi=300) for x in ("pdf","png")]
print(f"Tg: {TG}"); print(f"quiebre Tb = {Tb0} K (chi2+1 {min(IC)}-{max(IC)}), chi2 = {chi0:.1f}")
print("pico kappa(t) 2e10:",{T:(round(v[0],2),round(v[1],2),round(v[2]*1000)) for T,v in PICO.items()},"(W/m/K, sem, fs)")
print(f"900 K: eta {E900}, kappa {K900}, D {D900}; 800 K: eta(200 ps) {E800}, D {D800}")
print(f"Arrhenius 1000-1600 K: Ea = {p[0]*1000*8.617333e-5:.3f} eV; SE 900 K / liquido = {se[E2[:,0]==900][0]/ref:.3f}")
