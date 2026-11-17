"""Combina replicas: eta y kappa vs ventana, por replica y global (30 bloques).
Uso: venv_structural/bin/python analisis_replicas.py <T> <Vmean> <tag1> <tag2> ... """
import sys, numpy as np
T, V, tags = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3:]
kB = 1.380649e-23; eV = 1.602176634e-19; nblk = 10
def acf(x, nlag):
    n = len(x); f = np.fft.rfft(x, 2*n); return np.fft.irfft(f*np.conj(f))[:nlag]/(n-np.arange(nlag))
def blocks(cols, nlag):
    n = len(cols[0]); m = n//nblk
    return np.array([np.mean([acf(c[b*m:(b+1)*m], nlag) for c in cols], axis=0) for b in range(nblk)])
ETA, KAP = [], []
for tg in tags:
    P = np.loadtxt(f"{tg}_press.dat", comments="#"); dtp = (P[1,0]-P[0,0])*1e-3
    nl = int(20/dtp); B = blocks([P[:,4]*1e5, P[:,5]*1e5, P[:,6]*1e5], nl)
    ETA.append(V*1e-30/(kB*T)*1e-12*(np.cumsum(B,1)-0.5*B[:,:1])*dtp*1e3)
    J = np.loadtxt(f"{tg}_hflux.dat", comments="#"); dtj = (J[1,0]-J[0,0])*1e-3
    nj = int(20/dtj); BJ = blocks([J[:,i]*eV*1e-10/1e-12 for i in (1,2,3)], nj)
    KAP.append(3/(3*V*1e-30*kB*T**2)*1e-12*(np.cumsum(BJ,1)-0.5*BJ[:,:1])*dtj)
ETA, KAP = np.array(ETA), np.array(KAP)   # (rep, bloque, lag)
def rep(name, A, dt):
    print(f"# {name}: por replica (media de 10 bloques) y global (sem sobre {A.shape[0]*A.shape[1]} bloques)")
    for tw in [1, 2, 3, 5, 10]:
        k = int(tw/dt); r = A[:,:,k].mean(1); g = A[:,:,k].ravel()
        print(f"t={tw:2d} ps  reps: "+"  ".join(f"{x:.3f}" for x in r)+f"   global: {g.mean():.3f} +- {g.std(ddof=1)/np.sqrt(len(g)):.3f}   ({100*g.std(ddof=1)/np.sqrt(len(g))/g.mean():.1f}%)")
rep("eta (mPa s)", ETA, dtp); rep("kappa (W/m/K)", KAP, dtj)
# D por especie y E drift
for tg in tags:
    m = np.loadtxt(f"{tg}_msd.dat", comments="#"); t = m[:,0]*1e-3; h = len(t)//2
    D = [np.polyfit(t[h:], m[h:,i], 1)[0]/6*10 for i in (1,2,3)]   # pendiente A^2/ps /6 -> D; 1 A^2/ps = 1e-8 m2/s = 10 (1e-9 m2/s)
    th = np.loadtxt(f"{tg}_thermo.dat", comments="#")
    print(f"{tg}: D(Zr,Cu,Ag)[1e-9 m2/s]=", [round(d,3) for d in D], f" Tmean={th[:,1].mean():.2f}  Etot {th[0,4]:.1f}->{th[-1,4]:.1f}")
