"""Criterio de ruido: t_c = ultimo lag donde |C(t)| > N_SIG * sem(t) (sem entre 30 bloques) dentro de 0-TMAX ps.
Mas alla de t_c la ACF es indistinguible de 0 y solo agrega ruido; kappa = trapecio hasta t_c.
Se imprime tambien kappa(t) trapecio a varios t para ver la forma de la integral."""
import sys, numpy as np
kB=1.380649e-23; eV=1.602176634e-19; NB=10; TAGS=["cal_r1","cal_r2","cal_r3"]; TMAX=5.0
NS=float(sys.argv[2]) if len(sys.argv)>2 else 3.0; rate=sys.argv[1] if len(sys.argv)>1 else "2e10"
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
print(f"rate {rate}; criterio |C|>{NS} sem; kappa trapecio [W/m/K]")
print(" T(K) cruce0  t_c(ps)  kappa(t_c)       kappa a t=0.15 0.25 0.35 0.5 1 2 ps")
for T in (300,400,500,600,1000,1200):
    d=f"cal_{T}K_{rate}"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); BL=[]
    for tg in TAGS:
        J=np.loadtxt(f"{d}/{tg}_hflux.dat",comments="#"); dt=(J[1,0]-J[0,0])*1e-3; nl=int(TMAX/dt)
        X=[J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)]; m=len(X[0])//NB
        BL += [np.mean([acf(x[b*m:(b+1)*m],nl) for x in X],axis=0) for b in range(NB)]
    BL=np.array(BL); pref=dt/(V*1e-30*kB*T**2)*1e-12
    Cm=BL.mean(0); sem=BL.std(0,ddof=1)/np.sqrt(len(BL))
    sig=np.where(np.abs(Cm)>NS*sem)[0]; ic=int(sig.max()) if len(sig) else 0
    K=(np.cumsum(BL,axis=1)-0.5*BL[:,[0]])*pref; z=np.where(Cm<0)[0][0]*dt
    row=" ".join(f"{K[:,int(t/dt)].mean():.3f}" for t in (0.15,0.25,0.35,0.5,1,2))
    print(f"{T:5d} {z:6.3f}  {ic*dt:6.3f}  {K[:,ic].mean():.3f}+-{K[:,ic].std(ddof=1)/np.sqrt(len(K)):.3f}   {row}")
