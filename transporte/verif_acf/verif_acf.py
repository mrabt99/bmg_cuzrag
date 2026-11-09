"""Compara fix ave/correlate (verif_corr.dat) contra numpy sobre la misma serie (verif_series.dat).
Uso: /scratch/negishi/aaabt/venv_structural/bin/python verif_acf.py
"""
import numpy as np

S = np.loadtxt("verif_series.dat", comments="#")[:, 1:]      # pxy pxz pyz
n = len(S)

blocks, cur = [], None
for ln in open("verif_corr.dat"):
    if ln.startswith("#"): continue
    p = ln.split()
    if len(p) == 2:                                          # "timestep nventanas"
        cur = []; blocks.append(cur)
    else:
        cur.append([float(x) for x in p])
L = np.array(blocks[-1])                                     # ultimo bloque (ave running)
nl = len(L); cnt = L[:, 2]; C = L[:, 3:6]
print(f"muestras serie n={n}   bloques en verif_corr.dat={len(blocks)}   lags={nl}")

def acf(x, nlag):
    f = np.fft.rfft(x, 2 * n)
    return np.fft.irfft(f * np.conj(f))[:nlag]

k = np.arange(nl)
raw = np.array([acf(S[:, i], nl) for i in range(3)]).T        # suma_i x_i x_{i+k}
print("\nNcount de LAMMPS vs esperado:")
for kk in [0, 1, 10, 100, nl - 1]:
    print(f"  lag {kk:4d}: Ncount={cnt[kk]:.0f}   n-k={n-kk}   n={n}")

for name, norm in [("n-k (insesgada)", (n - k)[:, None]), ("n (sesgada)", n)]:
    ref = raw / norm
    rel = np.abs(C - ref) / np.abs(ref).clip(1e-30)
    print(f"\nLAMMPS vs numpy/{name}: max |dif rel| = {rel.max():.3e}   (lag 0: {rel[0].max():.3e}, lag {nl-1}: {rel[-1].max():.3e})")
print("\nlag0 LAMMPS:", C[0], "  numpy sum/n:", raw[0] / n)
