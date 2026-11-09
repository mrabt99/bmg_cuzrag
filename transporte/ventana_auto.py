"""Ventana de integracion objetiva para kappa (Green-Kubo): ventana automatica de Sokal (1989).
tau_int(M) = 1/2 + sum_{k=1..M} rho_k (rho = C/C0, en muestras); se elige el menor M con M >= c*tau_int(M) (c=6).
kappa = trapecio de C hasta M*dt. Error: sem entre 30 bloques (3 replicas x 10) a ESE M. Ademas se informa el
cruce por cero de la ACF y kappa a M, M/2 y 2M (sensibilidad a la eleccion).
Uso: python ventana_auto.py [rate]    (T = 300,400,500,600,1000,1200)"""
import sys, numpy as np
kB=1.380649e-23; eV=1.602176634e-19; NB=10; TAGS=["cal_r1","cal_r2","cal_r3"]; C_SOK=6.0
rate = sys.argv[1] if len(sys.argv)>1 else "2e10"
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
print(f"rate {rate}, c={C_SOK}.  kappa con regla del TRAPECIO (C(0) con medio peso)")
print(" T(K)  cruce0(ps)  M*(ps)  kappa(M*)        kappa(M*/2)  kappa(2M*)   [W/m/K]")
for T in (300,400,500,600,1000,1200):
    d=f"cal_{T}K_{rate}"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); BL=[]
    for tg in TAGS:
        J=np.loadtxt(f"{d}/{tg}_hflux.dat",comments="#"); dt=(J[1,0]-J[0,0])*1e-3; nl=int(30/dt)
        X=[J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)]; m=len(X[0])//NB
        BL += [np.mean([acf(x[b*m:(b+1)*m],nl) for x in X],axis=0) for b in range(NB)]
    BL=np.array(BL); pref=dt/(V*1e-30*kB*T**2)*1e-12
    cum = lambda B: (np.cumsum(B,axis=-1)-0.5*B[...,[0]])*pref      # trapecio acumulado
    Cm=BL.mean(0); rho=Cm/Cm[0]; tau=0.5+np.cumsum(rho[1:])          # tau_int(M), M=1..
    M=np.arange(1,len(tau)+1); ok=np.where(M>=C_SOK*tau)[0]; Ms=int(M[ok[0]]) if len(ok) else None
    z=np.where(Cm<0)[0]
    if Ms is None: print(f"{T:5d}  no converge"); continue
    K=cum(BL); f=lambda mm: (K[:,min(mm,nl-1)].mean(), K[:,min(mm,nl-1)].std(ddof=1)/np.sqrt(len(K)))
    a,b,c=f(Ms),f(max(Ms//2,1)),f(2*Ms)
    print(f"{T:5d}  {z[0]*dt:8.3f}  {Ms*dt:7.3f}  {a[0]:.3f}+-{a[1]:.3f}   {b[0]:.3f}+-{b[1]:.3f}   {c[0]:.3f}+-{c[1]:.3f}")
