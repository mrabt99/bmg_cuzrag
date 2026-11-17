"""Green-Kubo: viscosidad y conductividad termica desde las series de in.nve.
Uso: venv_structural/bin/python analisis_gk.py <tag> <T> <Vmean_A3> [tmax_ps] [nbloques]
Lee <tag>_press.dat (bar, cada sp fs) y <tag>_hflux.dat (eV*A/ps, SIN dividir por V, cada sj fs).
Escribe <tag>_gk.png y <tag>_gk.txt (integrales acumuladas) e imprime valores a distintas ventanas.
"""
import sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

tag, T, V = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
tmax = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0   # ps, ventana de correlacion
nblk = int(sys.argv[5]) if len(sys.argv) > 5 else 5
kB = 1.380649e-23; eV = 1.602176634e-19

def load(fn):
    return np.loadtxt(fn, comments="#")[:, 1:]          # columna 0 = paso

def acf(x, nlag):
    """autocorrelacion <x(0)x(t)> (sesgada-normalizada por n-k) via FFT, x: (n,) """
    n = len(x); f = np.fft.rfft(x - 0, 2 * n)
    c = np.fft.irfft(f * np.conj(f))[:nlag]
    return c / (n - np.arange(nlag))

def gk(series, dt_ps, nlag, nblk):
    """series: lista de arrays 1D (componentes). Devuelve ACF media y por bloque."""
    n = len(series[0]); m = n // nblk
    blk = []
    for b in range(nblk):
        blk.append(np.mean([acf(s[b*m:(b+1)*m], nlag) for s in series], axis=0))
    return np.array(blk)

# ---------- viscosidad ----------
P = load(f"{tag}_press.dat"); dtp = np.loadtxt(f"{tag}_press.dat", comments="#")[:2, 0]
dtp = (dtp[1] - dtp[0]) * 1e-3                               # ps (1 paso = 1 fs)
pxy, pxz, pyz = P[:, 3] * 1e5, P[:, 4] * 1e5, P[:, 5] * 1e5  # Pa
nl = int(tmax / dtp)
pref = V * 1e-30 / (kB * T) * 1e-12                          # m3/(J) * s por ps
res = {}
for name, comps in [("pxy", [pxy]), ("pxz", [pxz]), ("pyz", [pyz]), ("3 comp", [pxy, pxz, pyz])]:
    B = gk(comps, dtp, nl, nblk)
    I = pref * (np.cumsum(B, axis=1) - 0.5*B[:, :1]) * dtp                    # Pa s, acumulada por bloque
    res[name] = (B, I)
t = np.arange(nl) * dtp

# ---------- conductividad ----------
try:
    J = load(f"{tag}_hflux.dat"); dtj = np.loadtxt(f"{tag}_hflux.dat", comments="#")[:2, 0]
    dtj = (dtj[1] - dtj[0]) * 1e-3
    Jc = [J[:, i] * eV * 1e-10 / 1e-12 for i in range(3)]    # J*m/s (aun sin /V)
    nlj = int(tmax / dtj); BJ = gk(Jc, dtj, nlj, nblk)
    kpref = 1.0 / (3 * V * 1e-30 * kB * T**2) * 1e-12
    KJ = kpref * (np.cumsum(BJ, axis=1) - 0.5*BJ[:, :1]) * dtj * 3             # *3: acf(3 comps) ya promedia; sumar -> J.J
    tj = np.arange(nlj) * dtj
except OSError:
    KJ = None

# ---------- reporte ----------
print(f"# {tag}  T={T} K  V={V:.4g} A^3  dt_p={dtp} ps  ventana={tmax} ps  bloques={nblk}")
print("# integral acumulada de viscosidad (mPa s), media +- sem de bloques")
for tw in [0.1, 0.5, 1, 2, 5, 10, tmax - dtp]:
    if tw > t[-1]: continue
    k = int(tw / dtp)
    row = [f"t={tw:5.2f} ps"]
    for name, (B, I) in res.items():
        v = I[:, k] * 1e3; row.append(f"{name}: {v.mean():.3f}+-{v.std(ddof=1)/np.sqrt(nblk):.3f}")
    print("  ".join(row))
if KJ is not None:
    print("# integral acumulada de kappa (W/m/K)")
    for tw in [0.1, 0.5, 1, 2, 5, 10, tmax - dtj]:
        if tw > tj[-1]: continue
        k = int(tw / dtj); v = KJ[:, k]
        print(f"t={tw:5.2f} ps  kappa: {v.mean():.4f}+-{v.std(ddof=1)/np.sqrt(nblk):.4f}")

fig, ax = plt.subplots(2 if KJ is not None else 1, 2, figsize=(11, 7), squeeze=False)
for name, (B, I) in res.items():
    ax[0, 0].plot(t, B.mean(0) / B.mean(0)[0], label=name)
    ax[0, 1].plot(t, I.mean(0) * 1e3, label=name)
ax[0, 0].set(xscale="log", xlabel="t (ps)", ylabel="ACF normalizada", title="Presion"); ax[0, 0].axhline(0, c="k", lw=.5)
ax[0, 1].set(xscale="log", xlabel="t (ps)", ylabel="eta acumulada (mPa s)", title="Integral"); ax[0, 0].legend()
if KJ is not None:
    ax[1, 0].plot(tj, BJ.mean(0) / BJ.mean(0)[0]); ax[1, 0].set(xscale="log", xlabel="t (ps)", title="ACF flujo de calor")
    ax[1, 1].plot(tj, KJ.mean(0)); ax[1, 1].set(xscale="log", xlabel="t (ps)", ylabel="kappa (W/m/K)")
fig.tight_layout(); fig.savefig(f"{tag}_gk.png", dpi=120)
