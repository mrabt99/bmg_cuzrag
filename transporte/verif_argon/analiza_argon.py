"""Compara el kappa del ejemplo Ar de la doc de LAMMPS con nuestro acf/conversion.
real : serie en kcal/mol*A/fs; formula identica a la de la doc (scale con kCal2J, trap) -> debe dar lo que imprime LAMMPS.
metal: serie en eV*A/ps; conversion de nuestro pipeline (eV*1e-10/1e-12, kB=1.380649e-23, cumsum*dt, mean de 3 comps).
Uso: python analiza_argon.py   (en esta carpeta, tras terminar run_argon.sbatch)"""
import re, glob, numpy as np
T = 70.0; P = 200                      # lags (correlation length de la doc)
eV = 1.602176634e-19
def acf(x, nl):
    n = len(x); f = np.fft.rfft(x, 2*n); return np.fft.irfft(f*np.conj(f))[:nl]/(n-np.arange(nl))
def acf3(J, nl): return np.mean([acf(J[:, i], nl) for i in range(3)], axis=0)
def vol(log): return float(re.search(r"V_A3 = ([0-9.eE+-]+)", open(log).read()).group(1))
def lmp_k(log):
    m = re.search(r"average conductivity: ([0-9.eE+-]+)", open(log).read()); return float(m.group(1)) if m else np.nan
def trap(c, h): return h*(c.sum()-0.5*(c[0]+c[-1]))
res = {}
for ens in ("nvt", "nve"):
    for v in ("real", "metal"):
        ks = []
        for lg in sorted(glob.glob(f"log.ar_{v}_{ens}_*")):
            seed = lg.split("_")[-1]; V = vol(lg)
            D = np.loadtxt(f"ar_{v}_{ens}_{seed}_Jraw.dat", comments="#")
            J = D[:, 1:4]
            if v == "real":
                kB=1.3806504e-23; kc=4186.0/6.02214e23; conv=kc*kc/1e-15/1e-10; dt=4.0; s=10
                C = acf3(J, P+1)
                # la doc: trap sobre f_JJ (una componente, ave running) * scale; promedio de 3 comps = trap(mean C)
                k = trap(C, 1.0)*conv/kB/T/T/V*s*dt
                kl = lmp_k(lg); ks.append((seed, k, kl, np.nan))
            else:
                kB=1.380649e-23; h = 10*0.004      # ps entre muestras
                Js = J*eV*1e-10/1e-12
                C = acf3(Js, P+1)
                pref = 1.0/(V*1e-30*kB*T**2)*1e-12
                kc_ = np.cumsum(C)[-1]*h*pref      # como kappa_vidrio.py (rectangulo, C(0) completo)
                kt = trap(C, h)*pref
                ks.append((seed, kc_, kt, np.nan))
        res[(ens, v)] = ks
for (ens, v), ks in res.items():
    print(f"\n[{ens}] {v}")
    for s_, a, b, _ in ks:
        print(f"  seed {s_:>7}: " + (f"mi acf+doc = {a:.4f}   LAMMPS (trap, log) = {b:.4f}  W/m/K" if v=="real"
                                      else f"nuestro pipeline cumsum = {a:.4f}   trapecio = {b:.4f}  W/m/K"))
    a = np.array([x[1] for x in ks]); 
    if len(a)>1: print(f"  media {a.mean():.4f} +- {a.std(ddof=1)/np.sqrt(len(a)):.4f}  (sem entre {len(a)} semillas)")
