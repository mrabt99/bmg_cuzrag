"""[2026-11-04: integral por TRAPECIO (C(0) a medio peso), como el trap() de LAMMPS; antes era rectangulo, que sesga +6-18 %]
kappa por Green-Kubo en vidrio: cal_300K_<tasa>. Uso: venv_structural/bin/python kappa_vidrio.py [T] [tasa ...]"""
import sys, numpy as np
T = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
rates = sys.argv[2:] or ["2e10", "1e12"]
kB = 1.380649e-23; eV = 1.602176634e-19; NB = 10; TAGS = ["cal_r1", "cal_r2", "cal_r3"]; TMAX = 30.0
TW = [0.25, 0.5, 1, 2, 3, 5, 10, 20, 29]

def acf(x, nl):
    n = len(x); f = np.fft.rfft(x, 2*n)
    return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))

res = {}
for r in rates:
    d = f"cal_{int(T)}K_{r}"
    V = float(open(f"{d}/cal_vmean.inc").read().split()[3])
    K, Jm, Bmean = [], [], []
    for tg in TAGS:
        J = np.loadtxt(f"{d}/{tg}_hflux.dat", comments="#"); dt = (J[1,0]-J[0,0])*1e-3
        nl = int(TMAX/dt); X = [J[:, i]*eV*1e-10/1e-12 for i in (1,2,3)]
        Jm.append([x.mean()/x.std() for x in X])
        m = len(X[0])//NB
        B = np.array([np.mean([acf(x[b*m:(b+1)*m], nl) for x in X], axis=0) for b in range(NB)])
        K.append((np.cumsum(B,1)-0.5*B[:,:1])*dt/(V*1e-30*kB*T**2)*1e-12); Bmean.append(B.mean(0))
    K = np.concatenate(K); t = np.arange(nl)*dt
    Bm = np.mean(Bmean, axis=0); Bm /= Bm[0]
    th = [np.loadtxt(f"{d}/{tg}_thermo.dat", comments="#") for tg in TAGS]
    res[r] = (K, dt)
    print(f"\n######## {r}  V={V:.6g} A^3  dt_J={dt} ps")
    print("media/std de J por replica (Jx,Jy,Jz):", [[round(v,3) for v in a] for a in Jm])
    print("T media NVE por replica:", [round(x[:,1].mean(),2) for x in th], " Etot ini->fin:", [(round(x[0,4]), round(x[-1,4])) for x in th])
    print("Ep/atomo medio (eV):", [round(x[:,2].mean()/65536,5) for x in th])
    ze = np.where(Bm<0)[0]; e1 = np.where(Bm<np.exp(-1))[0]
    print(f"ACF de J: 1/e en {t[e1[0]]:.3f} ps; primer cruce por 0 en {t[ze[0]]:.3f} ps")
    print("kappa(t) [W/m/K], global (sem sobre 30 bloques) y por replica:")
    for tw in TW:
        k = int(tw/dt); g = K[:,k]; rep = [K[i*NB:(i+1)*NB,k].mean() for i in range(3)]
        print(f"  t={tw:5.2f} ps: {g.mean():.4f} +- {g.std(ddof=1)/np.sqrt(len(g)):.4f} ({100*g.std(ddof=1)/np.sqrt(len(g))/g.mean():.1f}%)   reps: "+"  ".join(f"{x:.3f}" for x in rep))

if len(rates) == 2:
    (Ka, dt), (Kb, _) = res[rates[0]], res[rates[1]]
    print(f"\n######## diferencia {rates[0]} - {rates[1]}")
    for tw in [1, 2, 3, 5, 10]:
        k = int(tw/dt); a, b = Ka[:,k], Kb[:,k]
        da = a.mean()-b.mean(); s = np.hypot(a.std(ddof=1)/np.sqrt(len(a)), b.std(ddof=1)/np.sqrt(len(b)))
        print(f"  t={tw:2d} ps: dif = {da:+.4f} +- {s:.4f}  ({da/s:+.1f} sigma, {100*da/b.mean():+.1f}% de {rates[1]})")
