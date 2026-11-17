"""VDOS por especie del vidrio 2e10 a 300 K a partir de la VACF (in.vdos, 2026-11-14).

Lee  : v300_2e10_vacf_<i>.dat (i = 1..NORIG), columnas: paso, Zr, Cu, Ag (total), Zr x/y/z
       v300_2e10_thermo.dat
Sale : vdos_300K_2e10.csv  (nu en THz; g_Zr, g_Cu, g_Ag, g_tot normalizadas a ∫g dnu = 1)
       vacf_300K_2e10.csv  (t en ps; VACF normalizada c(t) = C(t)/C(0) por especie, media y σ entre orígenes)
Método: g(nu) = 4 ∫_0^tc c(t) w(t) cos(2 pi nu t) dt, w = cos^2(pi t / 2 tc), tc = 3 ps
        (misma ventana que Re kappa(nu) en bitacora/figuras/fig_tr_kappa_frecuencia.py).
        g_tot = sum_a x_a g_a (cada átomo aporta 3 modos).
Chequeos: equipartición m<v^2>/3kB = T por especie; isotropía x/y/z; ∫g dnu = 1;
          nu_D efectiva de g ≈ 3 nu^2 / nu_D^3 a baja frecuencia (vs 4,4 THz de las elásticas).
Uso: /scratch/negishi/aaabt/venv_structural/bin/python analiza_vdos.py
"""
import glob
import re

import numpy as np

TAG = "v300_2e10"
T = 300.0
TC = 3.0
KB, AMU = 1.380649e-23, 1.66053906660e-27
MASA = {"Zr": 91.224, "Cu": 63.546, "Ag": 107.8682}
N = {"Zr": 31129, "Cu": 31130, "Ag": 3277}
ESP = ["Zr", "Cu", "Ag"]

arch = sorted(glob.glob(f"{TAG}_vacf_*.dat"), key=lambda s: int(re.findall(r"_(\d+)\.dat", s)[0]))
D = [np.loadtxt(a, comments="#") for a in arch]
nmin = min(len(d) for d in D)
D = np.array([d[:nmin] for d in D])            # (norig, nt, 8)
t = (D[0, :, 0] - D[0, 0, 0]) * 1e-3            # ps
print(f"chequeo: {len(arch)} orígenes, {nmin} puntos, dt = {(t[1]-t[0])*1e3:.0f} fs, t_max = {t[-1]:.2f} ps")

th = np.loadtxt(f"{TAG}_thermo.dat", comments="#")
print(f"chequeo: NVE T media = {th[:,1].mean():.1f} K (σ {th[:,1].std():.1f}); "
      f"deriva de E_tot = {(th[-1,3]-th[0,3]):.3f} eV en {(th[-1,0]-th[0,0])*1e-3:.0f} ps")

C = {e: D[:, :, 1 + k] for k, e in enumerate(ESP)}
for e in ESP:
    Te = MASA[e] * AMU * C[e][:, 0] * 1e4 / (3 * KB)
    print(f"chequeo: equipartición {e}: T = {Te.mean():.1f} ± {Te.std(ddof=1):.1f} K (entre orígenes)")
cx = D[:, 0, 4:7].mean(0) / D[:, 0, 1].mean()
print(f"chequeo: isotropía Zr <vx²>,<vy²>,<vz²> / <v²> = {np.round(cx, 4)} (esperado 1/3; la columna"
      " total de compute vacf es la suma)")

c = {e: C[e] / C[e][:, :1] for e in ESP}         # normalizada por origen
m = t <= TC
tt = t[m]
win = np.cos(np.pi * tt / (2 * TC)) ** 2
nu = np.linspace(0, 15, 1501)
ker = np.cos(2 * np.pi * nu[:, None] * tt[None, :])
g, gs = {}, {}
for e in ESP:
    G = np.array([4 * np.trapezoid(ci[m] * win * ker, tt, axis=1) for ci in c[e]])
    g[e], gs[e] = G.mean(0), G.std(0, ddof=1) / np.sqrt(len(G))
x = {e: N[e] / sum(N.values()) for e in ESP}
g["tot"] = sum(x[e] * g[e] for e in ESP)
gs["tot"] = np.sqrt(sum((x[e] * gs[e]) ** 2 for e in ESP))

nu_l = np.linspace(0, 100, 20001)
for e in ESP:
    gl = 4 * np.trapezoid(c[e].mean(0)[m] * win * np.cos(2 * np.pi * nu_l[:, None] * tt[None, :]), tt, axis=1)
    j = g[e].argmax()
    nmed = np.trapezoid(nu * g[e], nu) / np.trapezoid(g[e], nu)
    print(f"chequeo: {e}: ∫g = {np.trapezoid(gl, nu_l):.4f} (0-100 THz), ∫_0^15 = {np.trapezoid(g[e], nu):.4f}; "
          f"máximo a {nu[j]:.2f} THz; <nu> = {nmed:.2f} THz; g(0) = {g[e][0]:.4f} 1/THz")
j = g["tot"].argmax()
print(f"chequeo: total: máximo a {nu[j]:.2f} THz; g(0) = {g['tot'][0]:.4f} 1/THz "
      "(≠ 0 en un vidrio a 300 K = difusión/relajación residual + resolución de la ventana)")
for a, b in [(0.4, 0.8), (0.5, 1.0), (0.8, 1.2)]:
    s = (nu >= a) & (nu <= b)
    k = np.mean(g["tot"][s] / nu[s] ** 2)
    print(f"chequeo: g/nu² en {a}-{b} THz = {k:.4f} THz^-3 -> nu_D efectiva = {(3/k)**(1/3):.2f} THz "
          "(Debye de las elásticas: 4,4 THz; menor = exceso de modos blandos / pico bosónico o resolución)")

np.savetxt("vdos_300K_2e10.csv",
           np.column_stack([nu] + [g[e] for e in ESP + ["tot"]] + [gs[e] for e in ESP + ["tot"]]),
           delimiter=",", fmt="%.6g",
           header="nu_THz,g_Zr,g_Cu,g_Ag,g_tot,sem_Zr,sem_Cu,sem_Ag,sem_tot", comments="")
np.savetxt("vacf_300K_2e10.csv",
           np.column_stack([t] + [c[e].mean(0) for e in ESP] + [c[e].std(0, ddof=1) for e in ESP]),
           delimiter=",", fmt="%.6g",
           header="t_ps,c_Zr,c_Cu,c_Ag,sd_Zr,sd_Cu,sd_Ag", comments="")
print("escrito: vdos_300K_2e10.csv, vacf_300K_2e10.csv (sem = σ entre orígenes / √n)")
