"""Prueba del termino convectivo (README.md). 2e10 a 600, 900 y 1200 K, 3 replicas x 1 ns, bloques de 100 ps.
Flujos (por componente x,y,z; SI):
  J  = flujo de calor total (hflux.dat)        Jc = sum_i e_i v_i (conv.dat col 1-3)      Jv = J - Jc (virial)
  D_s = N_s v_cm,s (flujo de numero de la especie s; vcm.dat). Independientes: D_Zr y D_Ag (D_Cu fijado por momento nulo).
Matriz de Green-Kubo L_ab(t) = 1/(3 V kB) int_0^t <X_a(0).X_b(t')> dt' (correlacion cruzada simetrizada, trapecio, normalizacion n-k).
  kappa_total(t) = L_JJ/T^2                      particion: L_JcJc, L_JvJv, 2 L_JcJv
  kappa(c): J -> J + c*(D_Zr+D_Cu+D_Ag), corrimiento del cero de energia de todos los atomos en c eV (el kappa total depende de c)
  kappa_proy(t) = [L_JJ - L_JD L_DD^-1 L_DJ]/T^2 (flujo de calor a flujo de masa nulo; no depende de c)
Ventana 0.5-2 ps como kappa_tasas.py; curvas hasta 20 ps. Error: sem de 30 bloques (lineales), bootstrap sobre bloques (proyectado).
T: media de las replicas (difieren <= 3 K; efecto <= 1 %).
Uso (desde transporte/convectivo/): python analiza_conv.py -> analiza_conv.txt, analiza_conv.png"""
import re, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
kB=1.380649e-23; eV=1.602176634e-19; NB=10; TMAX=20.0; W0,W1=0.5,2.0; rng=np.random.default_rng(0)
M={"Zr":91.224,"Cu":63.546,"Ag":107.8682}; TS=(600,900,1200); CS=(-5,-2,-1,0,1,2,5)
SI=eV*1e-10/1e-12; SIn=1e-10/1e-12   # eV*A/ps -> J*m/s ; atomos*A/ps -> m/s
def xcorr(a,b,nl):
    n=len(a); fa=np.fft.rfft(a,2*n); fb=np.fft.rfft(b,2*n); c=np.fft.irfft(np.conj(fa)*fb)[:nl]/(n-np.arange(nl)); return c
def caso(T):
    d=f"conv_{T}K_2e10"; V=float(open(f"../cal_{T}K_2e10/cal_vmean.inc").read().split()[3])*1e-30; Ls=[]; Tm=[]; es=[]
    for r in (1,2,3):
        t=f"{d}/conv_r{r}"; log=open(glob.glob(f"{d}/log.conv_r{r}.*")[0]).read()
        N={g:int(re.search(rf"(\d+) atoms in group {g}",log).group(1)) for g in M}
        J=np.loadtxt(t+"_hflux.dat",comments="#"); C=np.loadtxt(t+"_conv.dat",comments="#"); Vc=np.loadtxt(t+"_vcm.dat",comments="#")
        Tm.append(np.loadtxt(t+"_thermo.dat",comments="#")[:,1].mean()); dt=(J[1,0]-J[0,0])*1e-3; nl=int(TMAX/dt)+1
        DZ=N["Zr"]*Vc[:,1:4]; DA=N["Ag"]*Vc[:,7:10]; DC=N["Cu"]*Vc[:,4:7]
        DC=-(M["Zr"]*DZ+M["Ag"]*DA)/M["Cu"]   # momento nulo exacto (vcm de salida redondeado)
        # energia media "efectiva" por especie: pendiente de la convectiva de la especie contra su flujo de numero
        es.append([np.sum(C[:,4+3*i:7+3*i]*Dd)/np.sum(Dd*Dd) for i,Dd in enumerate((DZ,N["Cu"]*Vc[:,4:7],DA))])
        X={"J":J[:,1:4]*SI,"Jc":C[:,1:4]*SI,"DZ":DZ*SIn,"DA":DA*SIn,"Dn":(DZ+DC+DA)*SIn}; X["Jv"]=X["J"]-X["Jc"]
        keys=list(X); m=len(J)//NB
        for b in range(NB):
            sl=slice(b*m,(b+1)*m); L={}
            for i,a in enumerate(keys):
                for c in keys[i:]:
                    cc=np.mean([0.5*(xcorr(X[a][sl,k],X[c][sl,k],nl)+xcorr(X[c][sl,k],X[a][sl,k],nl)) for k in range(3)],0)
                    L[a,c]=L[c,a]=(np.cumsum(cc)-0.5*cc[0])*dt*1e-12/(V*kB)
            Ls.append(L)
    return Ls,np.arange(nl)*dt,np.mean(Tm),np.mean(es,0)
def kproy(L,T):   # L: dict de arrays (por bloque o medio)
    LDD=np.array([[L["DZ","DZ"],L["DZ","DA"]],[L["DA","DZ"],L["DA","DA"]]]).transpose(2,0,1)
    LJD=np.array([L["J","DZ"],L["J","DA"]]).T
    sol=np.linalg.solve(LDD,LJD[...,None])[...,0]; return (L["J","J"]-np.sum(LJD*sol,1))/T**2
out=[]; R={}
for T in TS:
    Ls,t,Tm,es=caso(T); sel=(t>=W0)&(t<=W1); nb=len(Ls)
    mean=lambda k: np.mean([L[k] for L in Ls],0)
    def stat(f):   # f(L)->curva ; valor en ventana con sem de bloques
        g=np.array([f(L)[sel].mean() for L in Ls]); return g.mean(),g.std(ddof=1)/np.sqrt(nb)
    kt=stat(lambda L:L["J","J"]/Tm**2); kc=stat(lambda L:L["Jc","Jc"]/Tm**2); kv=stat(lambda L:L["Jv","Jv"]/Tm**2); kx=stat(lambda L:2*L["Jc","Jv"]/Tm**2)
    Lm={k:mean(k) for k in Ls[0]}; kp_curve=kproy(Lm,Tm)
    bs=[];
    for _ in range(300):
        idx=rng.integers(0,nb,nb); Lb={k:np.mean([Ls[i][k] for i in idx],0) for k in Ls[0]}; bs.append(kproy(Lb,Tm)[sel].mean())
    kp=(kp_curve[sel].mean(),np.std(bs))
    kc_=[stat(lambda L,c=c:(L["J","J"]+2*c*eV*L["J","Dn"]+(c*eV)**2*L["Dn","Dn"])/Tm**2) for c in CS]
    at=lambda f,x: f[int(round(x/(t[1]-t[0])))]
    R[T]=dict(t=t,tot=Lm["J","J"]/Tm**2,cc=Lm["Jc","Jc"]/Tm**2,vv=Lm["Jv","Jv"]/Tm**2,cv=2*Lm["Jc","Jv"]/Tm**2,proy=kp_curve,kc=kc_)
    out.append(f"===== {T} K (T media {Tm:.1f} K, {nb} bloques de 100 ps)")
    out.append(f"  energia efectiva por especie (pendiente de la convectiva contra el flujo de numero): Zr {es[0]:.3f}, Cu {es[1]:.3f}, Ag {es[2]:.3f} eV")
    out.append(f"  ventana [{W0},{W1}] ps:  total {kt[0]:.3f}+-{kt[1]:.3f} | virial-virial {kv[0]:.3f}+-{kv[1]:.3f} | conv-conv {kc[0]:.3f}+-{kc[1]:.3f} | cruzada {kx[0]:+.3f}+-{kx[1]:.3f} | PROYECTADO {kp[0]:.3f}+-{kp[1]:.3f}")
    out.append("  kappa(c) con el cero de energia corrido c eV: "+", ".join(f"c={c:+d}: {v:.3f}" for c,(v,e) in zip(CS,kc_)))
    out.append("  kappa(t) [total / proyectado / conv-conv / virial-virial] a 2, 5, 10, 20 ps: "
               +" | ".join(f"{x:g} ps: {at(R[T]['tot'],x):.3f}/{at(kp_curve,x):.3f}/{at(R[T]['cc'],x):.3f}/{at(R[T]['vv'],x):.3f}" for x in (2,5,10,19.9)))
out.append("")
out.append("referencia kappa_tasas.csv (r1-r3, otras semillas): 600 K 0.825+-0.019, 900 K 0.789+-0.018, 1200 K 0.879+-0.017")
open("analiza_conv.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
fig,ax=plt.subplots(1,4,figsize=(15,3.6))
for a,T in zip(ax[:3],TS):
    q=R[T]
    for k,lab,st in (("tot","total","k-"),("proy","proyectado (flujo de masa nulo)","C3-"),("vv","virial-virial","C0--"),("cc","conv-conv","C1--"),("cv","cruzada","C2:")):
        if k=="proy" and T==600: continue   # sin difusion L_DD ~ 0: proyeccion mal condicionada (vale el total)
        a.plot(q["t"],q[k],st,lw=1,label=lab)
    a.set_ylim(-0.05,1.5)
    a.axvspan(W0,W1,color="0.92"); a.set_xscale("log"); a.set_xlim(0.02,20); a.axhline(0,color="k",lw=0.4); a.set_title(f"{T} K"); a.set_xlabel("t (ps)"); a.set_ylabel(r"$\kappa(t)$ (W/m/K)")
ax[0].legend(fontsize=6)
for T in TS: ax[3].plot(CS,[v for v,e in R[T]["kc"]],"o-",label=f"{T} K")
ax[3].set_xlabel("corrimiento del cero de energia c (eV)"); ax[3].set_ylabel(r"$\kappa$ total, ventana 0.5-2 ps"); ax[3].legend(fontsize=7)
fig.tight_layout(); fig.savefig("analiza_conv.png",dpi=150)
