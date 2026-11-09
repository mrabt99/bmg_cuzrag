"""Aplica el protocolo del notebook de Seba (ML_Lammps_ZrCuAg_Structure_Paper_1.ipynb, celdas 9-11) a NUESTRAS series de J:
 - J muestreado cada 100 fs (Nevery=100) -> se submuestrea nuestra serie de 5 fs tomando 1 de cada 20
 - DACF de toda la serie (1 ns) por FFT, normalizada por n-k, promedio de 3 componentes, sin restar la media
 - kappa_running = prefactor * cumtrapz(DACF, t)
 - kappa = media de kappa_running entre el 10 % y el 80 % de los lags (i_start=0.1*N, i_end=0.8*N)
Tambien se calcula con la ventana 0.5-2 ps y paso 5 fs (nuestro metodo) y con el paso de 100 fs a tiempos cortos.
Uso (desde transporte/): python replica_metodo_seba.py [tasa]"""
import sys, numpy as np
kB=1.380649e-23; eV=1.602176634e-19; rate=sys.argv[1] if len(sys.argv)>1 else "2e10"
def acf(x):
    n=len(x); f=np.fft.rfft(x,2*n); return np.fft.irfft(f*np.conj(f))[:n]/(n-np.arange(n))
def cumtrap(y,dx): return np.concatenate([[0],np.cumsum(0.5*(y[1:]+y[:-1])*dx)])
print(f"rate {rate}. kappa [W/m/K] por replica y media (3 replicas); paper: 300:2.43 400:1.37 500:1.12 600:0.77 1000:0.94")
print(" T(K)  | metodo Seba (100 fs, media 10-80% de lags) por replica | media | nuestro (5 fs, 0.5-2 ps) | Seba-100fs pero ventana 0.3-2 ps")
for T in (300,400,500,600,1000):
    d=f"cal_{T}K_{rate}"; V=float(open(f"{d}/cal_vmean.inc").read().split()[3])*1e-30; ks=[]; kn=[]; kw=[]
    for tg in ("cal_r1","cal_r2","cal_r3"):
        J=np.loadtxt(f"{d}/{tg}_hflux.dat",comments="#"); dt5=(J[1,0]-J[0,0])*1e-3
        X=[J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)]
        # --- Seba: submuestreo a 100 fs ---
        s=int(round(0.1/dt5)); Xs=[x[::s] for x in X]; C=np.mean([acf(x) for x in Xs],axis=0)
        t=np.arange(len(C))*0.1e-12; kr=cumtrap(C,0.1e-12)/(V*kB*T**2)
        N=len(kr); ks.append(kr[int(0.1*N):int(0.8*N)].mean())
        # --- nuestro: 5 fs, ventana 0.5-2 ps, ACF sobre bloques de 100 ps simplificado: ACF global hasta 3 ps ---
        C5=np.mean([acf(x)[:int(3/dt5)] for x in X],axis=0); k5=cumtrap(C5,dt5*1e-12)/(V*kB*T**2)
        kn.append(k5[int(0.5/dt5):int(2/dt5)+1].mean())
        kw.append(kr[3:21].mean())   # lags 0.3-2.0 ps con paso 100 fs
    print(f"{T:5d}  | "+"  ".join(f"{k:.3f}" for k in ks)+f" | {np.mean(ks):.3f} +- {np.std(ks,ddof=1)/np.sqrt(3):.3f} | {np.mean(kn):.3f} | {np.mean(kw):.3f}")
