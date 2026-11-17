"""Compara J (heat/flux) con stress/atom vs centroid/stress/atom. Lee J_ambos.dat y dstress_max.dat.
Criterio: |dJ|max / rms(J) ~ 1e-12 o menos => identicos (redondeo); ~1e-3 o mas => diferencia real.
Uso: /scratch/negishi/aaabt/venv_structural/bin/python compara.py"""
import numpy as np
d = np.loadtxt("J_ambos.dat", comments="#")
J, Jc, K, Kc = d[:, 1:4], d[:, 4:7], d[:, 7:10], d[:, 10:13]
rms = np.sqrt((J ** 2).mean())
print(f"{len(d)} muestras (cada 5 pasos); rms(J) = {rms:.4g} eV*A/ps")
for nom, a, b in (("J total", J, Jc), ("J convectivo", K, Kc)):
    dif = np.abs(a - b)
    print(f"  {nom:13s}: max|dJ| = {dif.max():.3e}   max|dJ|/rms(J) = {dif.max() / rms:.3e}")
# integral de la ACF (sin normalizar) con las dos series: si J es igual, da igual
for nom, a in (("stress/atom", J), ("centroid", Jc)):
    c0 = (a * a).sum(axis=1).mean() / 3
    print(f"  <J^2>/3 con {nom:11s} = {c0:.6e}")
s = np.loadtxt("dstress_max.dat", comments="#")
s = np.atleast_2d(s)
print("Tensor por atomo (bar*A^3), maximo sobre atomos y sobre los instantes guardados:")
print(f"  max|S_st - S_cst| (xx yy zz xy xz yz) = " + " ".join(f"{x:.2e}" for x in s[:, 1:7].max(axis=0)))
print(f"  asimetria centroide max|xy - yx|     = {s[:, 7].max():.2e}")
print(f"  escala: max(|Sxx|+|Syy|+|Szz|)        = {s[:, 8].max():.2e}")
