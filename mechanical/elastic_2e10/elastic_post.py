#!/usr/bin/env python3
"""Post-proceso de las constantes elasticas (300 K o 0 K). Lee Ccol1..6.txt (columnas de Cij en GPa; salida de in.elastic_dir) y, a 300 K, las series
ser_dir<j>_s<signo>.dat (promedios de presion cada 1 ps) para estimar el error; calcula:
  - matriz Cij 6x6 (filas = tension, columnas = deformacion, orden de Voigt xx yy zz yz xz xy), asimetria Cij vs Cji, simetrizada
  - isotropia: C11-C12-2*C44 (promediados), indice de anisotropia universal A_U
  - promedios isotropos de Voigt, Reuss y Hill: K, G, E, nu
  - velocidades v_L = sqrt((K + 4G/3)/rho), v_T = sqrt(G/rho), v_m = [(1/3)(1/v_L^3 + 2/v_T^3)]^(-1/3) y temperatura de Debye de Anderson
        Theta_D = (hbar/kB) (6 pi^2 n)^(1/3) v_m,   n = atomos/volumen   (comparar con los 312 K de Mendez)
Errores: a 300 K, error estandar por bloques de 5 ps de cada presion media -> propagacion Monte Carlo (1000 muestras) sobre Cij -> K, G, E, nu, Theta_D.
SOLO ruido estadistico de UNA corrida (no incluye variacion entre vidrios ni incertidumbre del potencial).
Uso: python elastic_post.py --dir 300K --up 0.005 --vol <A^3 o 'auto'>   |   python elastic_post.py --dir 0K/up_1e-4 --up 1e-4 --vol-from 0K/prep/log.prep
Salidas: <dir>/elastic_resumen.txt, <dir>/Cij.csv"""
import sys, re, os, glob
import numpy as np
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
D = arg("--dir"); UP = float(arg("--up")); CF = 1e-4                       # bar -> GPa
NAT = {1: 31129, 2: 31130, 3: 3277}; MASS = {1: 91.22, 2: 63.55, 3: 107.87}   # tipos 1 Zr, 2 Cu, 3 Ag (N_Zr = 31129, N_Cu = 31130, N_Ag = 3277)
AMU, HB, KB, NA = 1.66053907e-27, 1.054571817e-34, 1.380649e-23, 6.02214076e23

def read_cols(d):
    C = np.full((6, 6), np.nan)
    for j in range(1, 7):
        t = open(f"{d}/Ccol{j}.txt").read()
        for i, jj, v in re.findall(r"C(\d)(\d) = (-?[\d.]+)", t):
            assert int(jj) == j; C[int(i) - 1, j - 1] = float(v)
    return C

def read_vol(src):                                   # volumen de referencia (A^3) del log de preparacion
    t = open(src).read()
    m = re.findall(r"vol = ([\d.]+) A\^3", t)
    if m: return float(m[-1])                        # 300 K: 'caja fijada al promedio ... vol = X A^3'
    m = re.findall(r"caja Lx = ([\d.]+) Ly = ([\d.]+) Lz = ([\d.]+)", t)   # 0 K: caja de la referencia
    return float(np.prod([float(x) for x in m[-1]]))

def series_sem(d, j):
    """Error estandar (GPa) de cada Cij de la columna j a partir de los bloques de 5 ps de las series de presion."""
    out = np.zeros(6); order = [1, 2, 3, 6, 5, 4]    # columnas del .dat (pxx pyy pzz pxy pxz pyz) -> orden de Voigt xx yy zz yz xz xy
    sem = {}
    for s in (-1, 1):
        a = np.loadtxt(f"{d}/ser_dir{j}_s{s}.dat", comments="#")[:, 1:]; a = a[:, [o - 1 for o in order]]
        nb = max(a.shape[0] // 5, 1); b = a[:nb * 5].reshape(nb, -1, 6).mean(1) if a.shape[0] >= 10 else a
        sem[s] = b.std(0, ddof=1) / np.sqrt(len(b))
    return np.hypot(sem[-1], sem[1]) / (2 * UP) * CF

def moduli(C):
    Cs = 0.5 * (C + C.T); S = np.linalg.inv(Cs)
    d = np.diag(Cs)[:3].sum(); o = (Cs[0, 1] + Cs[0, 2] + Cs[1, 2]); sh = np.diag(Cs)[3:].sum()
    KV = (d + 2 * o) / 9; GV = (d - o + 3 * sh) / 15
    sd = np.diag(S)[:3].sum(); so = S[0, 1] + S[0, 2] + S[1, 2]; ssh = np.diag(S)[3:].sum()
    KR = 1 / (sd + 2 * so); GR = 15 / (4 * sd - 4 * so + 3 * ssh)
    return KV, GV, KR, GR

def derived(K, G, rho, n):
    E = 9 * K * G / (3 * K + G); nu = (3 * K - 2 * G) / (2 * (3 * K + G))
    vL = np.sqrt((K + 4 * G / 3) * 1e9 / rho); vT = np.sqrt(G * 1e9 / rho)
    vm = ((1 / 3) * (1 / vL ** 3 + 2 / vT ** 3)) ** (-1 / 3); th = HB / KB * (6 * np.pi ** 2 * n) ** (1 / 3) * vm
    return E, nu, vL, vT, vm, th

C = read_cols(D); vol = float(arg("--vol")) if arg("--vol") not in (None, "auto") else read_vol(arg("--vol-from", f"{D}/log.prep"))
Natoms = sum(NAT.values()); mass = sum(NAT[t] * MASS[t] for t in NAT) * AMU          # kg
rho = mass / (vol * 1e-30); n = Natoms / (vol * 1e-30)
Cs = 0.5 * (C + C.T); asym = np.abs(C - C.T).max()
has_ser = all(os.path.exists(f"{D}/ser_dir{j}_s{s}.dat") for j in range(1, 7) for s in (-1, 1))
E_ = np.zeros((6, 6))
if has_ser:
    for j in range(1, 7): E_[:, j - 1] = series_sem(D, j)
KV, GV, KR, GR = moduli(C); K, G = 0.5 * (KV + KR), 0.5 * (GV + GR)
lines = [f"dir = {D} ; up = {UP} ; volumen = {vol:.1f} A^3 ; rho = {rho/1e3:.4f} g/cc ; n = {n:.4e} atomos/m^3",
         "Cij [GPa] (filas: tension xx yy zz yz xz xy; columnas: deformacion), tal como salen de las columnas:"]
lines += ["   " + " ".join(f"{x:9.2f}" for x in row) for row in C]
if has_ser: lines += ["error estandar (bloques de 5 ps) [GPa]:"] + ["   " + " ".join(f"{x:9.2f}" for x in row) for row in E_]
lines += [f"asimetria maxima |Cij - Cji| = {asym:.2f} GPa (la matriz se simetriza)",
          f"isotropia: <C11,C22,C33> = {np.diag(Cs)[:3].mean():.2f} ; <C12,C13,C23> = {(Cs[0,1]+Cs[0,2]+Cs[1,2])/3:.2f} ; <C44,C55,C66> = {np.diag(Cs)[3:].mean():.2f} ; "
          f"(C11-C12)/2 = {(np.diag(Cs)[:3].mean() - (Cs[0,1]+Cs[0,2]+Cs[1,2])/3)/2:.2f} (debe coincidir con C44 si es isotropo) ; A_U = {5*GV/GR + KV/KR - 6:.4f}"]
rows = {}
for nm, (Kx, Gx) in {"Voigt": (KV, GV), "Reuss": (KR, GR), "Hill": (K, G)}.items():
    E, nu, vL, vT, vm, th = derived(Kx, Gx, rho, n); rows[nm] = (Kx, Gx, E, nu, vL, vT, vm, th)
    lines.append(f"{nm:5s}: K = {Kx:7.2f}  G = {Gx:6.2f}  E = {E:6.2f} GPa  nu = {nu:.4f}  vL = {vL:6.0f}  vT = {vT:6.0f}  vm = {vm:6.0f} m/s  Theta_D = {th:6.1f} K")
if has_ser:
    rng = np.random.default_rng(1); sm = []
    for _ in range(1000):
        Cr = C + rng.normal(size=(6, 6)) * E_; Kr, Gr = np.array(moduli(Cr)).reshape(2, 2).mean(0)
        sm.append((Kr, Gr) + derived(Kr, Gr, rho, n)[:2] + (derived(Kr, Gr, rho, n)[5],))
    sm = np.array(sm); lines.append("Hill, error estandar (Monte Carlo sobre los errores de Cij): " + "  ".join(f"{nm} = {sm[:, k].std():.2f}" for k, nm in enumerate(["K", "G", "E", "nu", "Theta_D"])))
lines.append("Referencia: Theta_D de Mendez = 312 K ; piloto de traccion 2e10 (300 K, 1e9 1/s): E = 54.6 GPa, nu ~ 0.41")
print("\n".join(lines)); open(f"{D}/elastic_resumen.txt", "w").write("\n".join(lines) + "\n")
np.savetxt(f"{D}/Cij.csv", C, delimiter=",", fmt="%.4f", header="Cij GPa; filas tension (xx yy zz yz xz xy), columnas deformacion")
