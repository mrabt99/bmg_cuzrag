"""Ajuste VFT eta = eta0*exp(D*T0/(T-T0)) a eta(T) a 2e10 (eta_T.csv + 900 K de eta_largo.py), ponderado por 1/sem.
Compara con Arrhenius y con la cota de 800 K (>= 300 mPa s). Fragilidad m = dlog10(eta)/d(Tg/T) en Tg.
Uso (desde transporte/): python eta_vft.py -> eta_vft.txt, eta_vft.png"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.genfromtxt("eta_T.csv",delimiter=",",names=True)
T=np.append(d["T_K"],900.0); eta=np.append(d["eta_mPas"],42.3); e=np.append(d["sem"],2.7)
o=np.argsort(T); T,eta,e=T[o],eta[o],e[o]; y=np.log(eta); w=eta/e   # sigma de ln(eta) = e/eta
best=None
for T0 in np.arange(0,880,1.0):          # barrido en T0, ajuste lineal de ln(eta) vs 1/(T-T0)
    x=1/(T-T0); p=np.polyfit(x,y,1,w=w); chi=np.sum((w*(y-np.polyval(p,x)))**2)
    if best is None or chi<best[0]: best=(chi,T0,p)
chi,T0,p=best; B=p[0]; eta0=np.exp(p[1]); dof=len(T)-3
xa=1/T; pa=np.polyfit(xa,y,1,w=w); chia=np.sum((w*(y-np.polyval(pa,xa)))**2)
# incertidumbre de T0: rango con chi2 <= chi2_min + 1 (re-escalado por chi2/dof si > 1)
s=max(chi/dof,1.0); ok=[]
for t0 in np.arange(0,880,1.0):
    x=1/(T-t0); q=np.polyfit(x,y,1,w=w)
    if np.sum((w*(y-np.polyval(q,x)))**2)<=chi+s: ok.append(t0)
f=lambda TT: eta0*np.exp(B/(TT-T0))
Tg=700.0   # Tg aproximada de nuestro vidrio a 2e10 (ver analisis Tg), solo para m
eta_g=f(Tg); m=B*Tg/(Tg-T0)**2/np.log(10)
out=[f"VFT: eta0 = {eta0:.3f} mPa s, B = D*T0 = {B:.0f} K, T0 = {T0:.0f} K (rango chi2+1: {min(ok):.0f}-{max(ok):.0f} K), D = {B/T0 if T0>0 else np.inf:.1f}",
     f"chi2/dof VFT = {chi:.1f}/{dof}; Arrhenius (Ea = {pa[0]*8.617333e-5:.3f} eV) chi2/dof = {chia:.1f}/{len(T)-2}",
     f"predicciones VFT: 800 K -> {f(800):.0f} mPa s (cota medida >= 300), 700 K -> {f(700):.2e} mPa s",
     f"fragilidad m en Tg={Tg:.0f} K (solo orientativa, extrapolada): {m:.0f}",
     "residuos (eta_VFT/eta - 1): "+", ".join(f"{t:.0f}:{f(t)/v-1:+.3f}" for t,v in zip(T,eta))]
fig,ax=plt.subplots(figsize=(5.2,3.8)); TT=np.linspace(780,1650,200)
ax.errorbar(1000/T,eta,e,fmt="o",label="GK 2e10"); ax.plot(1000/TT,f(TT),label=f"VFT, T0={T0:.0f} K")
ax.plot(1000/TT,np.exp(np.polyval(pa,1/TT)),"--",label="Arrhenius")
ax.annotate("",xy=(1000/800,600),xytext=(1000/800,300),arrowprops=dict(arrowstyle="->")); ax.plot([1000/800],[300],"v",color="k",label="800 K: cota inferior")
ax.set_yscale("log"); ax.set_xlabel("1000/T (1/K)"); ax.set_ylabel(r"$\eta$ (mPa s)"); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig("eta_vft.png",dpi=150)
open("eta_vft.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
