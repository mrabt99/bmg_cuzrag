"""Regla tipo Zhang-Otani-Maginn: corte t_c = primer t (>= t_min) donde std/|media| del acumulado entre los 30 bloques de 100 ps llega a FR (30,40,50 %).
Valor = promedio del acumulado medio entre t_c/2 y t_c. Trapecio. Uso: python ventana_zhang.py   (kappa en 300-1200 K, 2e10; eta en 1200 K, 2e10)"""
import numpy as np
kB=1.380649e-23; eV=1.602176634e-19; NB=10; TAGS=["cal_r1","cal_r2","cal_r3"]
def acf(x,nl):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
def cum(B): return np.cumsum(B,1)-0.5*B[:,:1]
def run(name,K,dt,tmin):
    t=np.arange(K.shape[1])*dt; m=K.mean(0); r=K.std(0,ddof=1)/np.abs(m); out=[]
    for fr in (0.3,0.4,0.5):
        i=np.where((t>=tmin)&(r>=fr))[0]; ic=int(i[0]) if len(i) else len(t)-1
        sel=(t>=t[ic]/2)&(t<=t[ic]); out.append(f"{100*fr:.0f}%: t_c={t[ic]:.2f} ps, valor[tc/2,tc]={K[:,sel].mean():.3f}")
    print(f"{name:12s} "+" | ".join(out))
print("kappa W/m/K (tmin=0.3 ps, despues del lobulo negativo)")
for T in (300,400,500,600,1000,1200):
    d=f"cal_{T}K_2e10"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); K=[]
    for tg in TAGS:
        J=np.loadtxt(f"{d}/{tg}_hflux.dat",comments="#"); dt=(J[1,0]-J[0,0])*1e-3; nl=int(30/dt)
        X=[J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)]; m=len(X[0])//NB
        B=np.array([np.mean([acf(x[b*m:(b+1)*m],nl) for x in X],axis=0) for b in range(NB)])
        K.append(cum(B)*dt/(V*1e-30*kB*T**2)*1e-12)
    run(f"{T} K",np.concatenate(K),dt,0.3)
print("eta mPa s, 1200 K (tmin=1 ps)")
d="cal_1200K_2e10"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3]); E=[]; T=1200.
for tg in TAGS:
    P=np.loadtxt(f"{d}/{tg}_press.dat",comments="#"); dt=(P[1,0]-P[0,0])*1e-3; nl=int(30/dt); m=(len(P))//NB
    B=np.array([np.mean([acf(P[b*m:(b+1)*m,c]*1e5,nl) for c in (4,5,6)],axis=0) for b in range(NB)])
    E.append(cum(B)*V*1e-30/(kB*T)*1e-12*dt*1e3)
run("1200 K eta",np.concatenate(E),dt,1.0)
