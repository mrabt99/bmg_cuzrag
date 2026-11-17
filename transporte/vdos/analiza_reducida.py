"""Densidad de estados reducida g(nu)/nu^2 respecto del nivel de Debye elástico (búsqueda del pico bosónico).

El pico bosónico se define como un MÁXIMO de g(nu)/nu^2 por encima del nivel de Debye 3/nu_D^3.
Un exceso sin máximo no alcanza para llamarlo así.

Lee  : <tag>_vacf_<i>.dat (salida de in.vdos; columnas: paso, Zr, Cu, Ag (total), Zr x/y/z)
Sale : reducida_<tag>.csv  (nu en THz; g_tot, g/nu^2 / Debye para cada tc; sem entre orígenes)
Método: c_tot(t) = sum_a x_a C_a(t)/C_a(0); g(nu) = 4 ∫_0^tc c(t) w(t) cos(2 pi nu t) dt,
        w = cos^2(pi t / 2 tc). Resolución (FWHM de la ventana) ≈ 0,72 / tc.
Debye: nu_D = 4,36 THz (Theta_D = 209,3 K, mechanical/elastic_2e10/300K/elastic_resumen.txt,
       v_m = 1827 m/s, n = 5,7018e28 m^-3; mismo vidrio y mismo volumen).
Chequeos: el máximo (si lo hay) tiene que estar por encima de la resolución y no moverse con tc;
          altura y frecuencia con su sem entre orígenes; con el tag de 4 ps reproduce la tabla del
          README (2,2 / 2,05 / 1,85 / 1,3 a 0,5 / 0,7 / 0,9 / 1,5 THz).
Uso: venv_structural/bin/python analiza_reducida.py <tag> [tc1 tc2 ...]
     p. ej. analiza_reducida.py v300_2e10 1.5 2 3 3.9 ; analiza_reducida.py v300_2e10_L20 3 5 10 15 19
"""
import glob
import re
import sys

import numpy as np

TAG = sys.argv[1] if len(sys.argv) > 1 else "v300_2e10"
TCS = [float(v) for v in sys.argv[2:]] or [3.0]
NU_D = 4.36
DEB = 3 / NU_D**3
N = np.array([31129, 31130, 3277])
X = N / N.sum()

arch = sorted(glob.glob(f"{TAG}_vacf_*.dat"), key=lambda s: int(re.findall(r"_(\d+)\.dat$", s)[0]))
D = [np.loadtxt(a, comments="#") for a in arch]
nmin = min(len(d) for d in D)
D = np.array([d[:nmin] for d in D])
t = (D[0, :, 0] - D[0, 0, 0]) * 1e-3
c = (D[:, :, 1:4] / D[:, :1, 1:4] * X).sum(2)        # (norig, nt) c_tot normalizada por origen
print(f"chequeo: {TAG}: {len(arch)} orígenes, t_max = {t[-1]:.2f} ps, dt = {(t[1]-t[0])*1e3:.0f} fs; "
      f"nivel de Debye 3/nu_D^3 = {DEB:.4f} THz^-3 (nu_D = {NU_D} THz)")
cola = c[:, t > min(5.0, 0.5 * t[-1])]
print(f"chequeo: cola de c(t) (t > {min(5.0, 0.5*t[-1]):.1f} ps): media {cola.mean():.2e}, sd {cola.std():.2e}")

nu = np.arange(0.02, 3.001, 0.01)
fil = [nu]
cab = ["nu_THz"]
for tc in TCS:
    if tc > t[-1] + 1e-9:
        print(f"aviso: tc = {tc} ps > t_max, se omite")
        continue
    m = t <= tc
    tt = t[m]
    w = np.cos(np.pi * tt / (2 * tc)) ** 2
    ker = np.cos(2 * np.pi * nu[:, None] * tt[None, :])
    G = np.array([4 * np.trapezoid(ci[m] * w * ker, tt, axis=1) for ci in c])   # por origen
    R = G / nu**2 / DEB
    Rm, Rs = R.mean(0), R.std(0, ddof=1) / np.sqrt(len(R))
    res = 0.72 / tc
    ok = nu >= 2 * res
    # máximo interior: punto con vecinos menores, por encima de 2x resolución
    loc = np.where(ok[1:-1] & (Rm[1:-1] > Rm[:-2]) & (Rm[1:-1] > Rm[2:]))[0] + 1
    txt = ", ".join(f"{nu[i]:.2f} THz (altura {Rm[i]:.2f} ± {Rs[i]:.2f})" for i in loc[:3]) or "ninguno"
    vals = " ".join(f"{Rm[np.argmin(abs(nu-v))]:.2f}±{Rs[np.argmin(abs(nu-v))]:.2f}"
                    for v in (0.2, 0.3, 0.5, 0.7, 0.9, 1.5, 2.0))
    print(f"tc = {tc:5.1f} ps (res ≈ {res:.3f} THz): g/nu²/Debye a 0,2/0,3/0,5/0,7/0,9/1,5/2,0 THz = {vals}")
    print(f"          máximos locales por encima de {2*res:.2f} THz: {txt}")
    fil += [G.mean(0), Rm, Rs]
    cab += [f"g_tc{tc:g}", f"red_tc{tc:g}", f"sem_red_tc{tc:g}"]

np.savetxt(f"reducida_{TAG}.csv", np.column_stack(fil), delimiter=",", fmt="%.6g",
           header=",".join(cab), comments="")
print(f"escrito: reducida_{TAG}.csv (red = g/nu²/(3/nu_D³); sem = σ entre orígenes / √n)")
