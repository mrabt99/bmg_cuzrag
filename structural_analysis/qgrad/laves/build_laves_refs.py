#!/usr/bin/env python3
"""Construye Cu2Zr Laves ideal como data files de LAMMPS (type 1=Zr, 2=Cu; el 3=Ag se declara sin atomos: mismo orden que
pair_coeff * * ZrCuAg.lammps.eam Zr Cu Ag). Solo numpy.
  C15 (MgCu2, Fd-3m, 24 at/celda): Zr(A) en 8a = fcc + (0,0,0),(1/4,1/4,1/4)  [subred diamante cubica]
                                  Cu(B) en 16d = fcc + (5/8,5/8,5/8),(3/8,7/8,1/8),(7/8,1/8,3/8),(1/8,3/8,7/8)
  C14 (MgZn2, P6_3/mmc, 12 at/celda), parametros ideales c/a=sqrt(8/3), z=1/16, x=5/6:
        Zr(A) en 4f: (1/3,2/3,z),(2/3,1/3,z+1/2),(2/3,1/3,-z),(1/3,2/3,-z+1/2)   [subred diamante hexagonal]
        Cu(B) en 2a: (0,0,0),(0,0,1/2)   y 6h: (x,2x,1/4),(-2x,-x,1/4),(x,-x,1/4),(-x,-2x,3/4),(2x,x,3/4),(-x,x,3/4)
Comprueba coordinacion por distancia: en Laves ideal Zr tiene 16 vecinos (12 Cu + 4 Zr) y Cu tiene 12 (6 Cu + 6 Zr).
Parametros de red iniciales (se relajan despues con el EAM real)."""
import numpy as np
from itertools import product

def write_data(path, cell, frac, types, comment):
    """cell: 3x3 (filas = vectores a1,a2,a3). Caja triclinica de LAMMPS via a,b,c (a1 || x, a2 en xy)."""
    a1, a2, a3 = map(np.array, cell)
    ax = np.linalg.norm(a1); bx = a2 @ a1 / ax; by = np.sqrt(a2 @ a2 - bx ** 2)
    cx = a3 @ a1 / ax; cy = (a2 @ a3 - bx * cx) / by; cz = np.sqrt(a3 @ a3 - cx ** 2 - cy ** 2)
    H = np.array([[ax, bx, cx], [0, by, cy], [0, 0, cz]])            # columnas = vectores
    xyz = (H @ np.array(frac).T).T
    with open(path, "w") as f:
        f.write(f"{comment}\n\n{len(frac)} atoms\n3 atom types\n\n0 {ax:.8f} xlo xhi\n0 {by:.8f} ylo yhi\n0 {cz:.8f} zlo zhi\n")
        if abs(bx) > 1e-9 or abs(cx) > 1e-9 or abs(cy) > 1e-9: f.write(f"{bx:.8f} {cx:.8f} {cy:.8f} xy xz yz\n")
        f.write("\nMasses\n\n1 91.224\n2 63.546\n3 107.868\n\nAtoms # atomic\n\n")
        for i, (t, p) in enumerate(zip(types, xyz), 1): f.write(f"{i} {t} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n")
    return H

def wrap_unique(fr, ty):
    seen, F, T = {}, [], []
    for p, t in zip(fr, ty):
        q = tuple(np.round(np.mod(p, 1.0), 6) % 1.0); q = tuple(0.0 if abs(x - 1) < 1e-6 else x for x in q)
        if q in seen: assert seen[q] == t; continue
        seen[q] = t; F.append(q); T.append(t)
    return np.array(F), T

fcc = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
# ---- C15
A = [(0, 0, 0), (.25, .25, .25)]; B = [(5/8, 5/8, 5/8), (3/8, 7/8, 1/8), (7/8, 1/8, 3/8), (1/8, 3/8, 7/8)]
fr, ty = [], []
for t, sites in ((1, A), (2, B)):
    for s in sites:
        for tr in fcc: fr.append(np.array(s) + tr); ty.append(t)
fr, ty = wrap_unique(fr, ty); assert len(fr) == 24 and ty.count(1) == 8 and ty.count(2) == 16
a15 = 7.40
NR15 = (3, 3, 3)
def supercell(fr, ty, n):
    F, T = [], []
    for i, j, k in product(range(n[0]), range(n[1]), range(n[2])):
        for p, t in zip(fr, ty): F.append((p + np.array([i, j, k])) / np.array(n)); T.append(t)
    return np.array(F), T
F, T = supercell(fr, ty, NR15)
H = write_data("refs/C15_init.data", np.diag([a15 * n for n in NR15]), F, T, "C15 Cu2Zr (MgCu2) 3x3x3, a0=7.40")
# ---- C14
a14 = 5.20; c14 = a14 * np.sqrt(8 / 3); z, x = 1 / 16, 5 / 6
Zr = [(1/3, 2/3, z), (2/3, 1/3, z + .5), (2/3, 1/3, -z), (1/3, 2/3, -z + .5)]
Cu = [(0, 0, 0), (0, 0, .5), (x, 2*x, .25), (-2*x, -x, .25), (x, -x, .25), (-x, -2*x, .75), (2*x, x, .75), (-x, x, .75)]
fr2 = [np.array(p) for p in Zr + Cu]; ty2 = [1] * 4 + [2] * 8
fr2, ty2 = wrap_unique(fr2, ty2); assert len(fr2) == 12 and ty2.count(1) == 4
cell14 = [(a14, 0, 0), (-a14 / 2, a14 * np.sqrt(3) / 2, 0), (0, 0, c14)]
n14 = (5, 3, 3)
def supercell_gen(fr, ty, n):
    F, T = [], []
    for i, j, k in product(range(n[0]), range(n[1]), range(n[2])):
        for p, t in zip(fr, ty): F.append((p + np.array([i, j, k])) / np.array(n)); T.append(t)
    return np.array(F), T
F2, T2 = supercell_gen(fr2, ty2, n14)
cell14s = [tuple(np.array(v) * n) for v, n in zip(cell14, n14)]
write_data("refs/C14_init.data", cell14s, F2, T2, "C14 Cu2Zr (MgZn2) 5x3x3, a0=5.20 c/a ideal")

# ---- coordinacion por distancia, enumerando imagenes periodicas (-2..2)^3 (valido tambien con celdas chicas / hexagonal)
def coord_check(cell, frac, types, rc):
    cell = np.array(cell, float); frac = np.array(frac); types = list(types); n = len(frac)
    shifts = np.array(list(product(range(-2, 3), repeat=3)))
    out = []
    for i in range(n):
        d = (frac[None, :, :] + shifts[:, None, :] - frac[i][None, None, :]).reshape(-1, 3) @ cell
        dist = np.linalg.norm(d, axis=1); ids = np.tile(np.arange(n), len(shifts))
        m = (dist > 0.1) & (dist < rc)
        out.append((int(m.sum()), int(sum(types[j] == 1 for j in ids[m]))))
    return out
def report(name, cell, fr_, ty_, rc):
    c = coord_check(cell, fr_, ty_, rc)
    zc = sorted({c[i] for i in range(len(c)) if ty_[i] == 1}); cc = sorted({c[i] for i in range(len(c)) if ty_[i] == 2})
    print(f"{name} (rc={rc} A): Zr (N_vecinos, N_Zr) = {zc} | Cu = {cc}   esperado Zr (16,4), Cu (12,6)")
    return zc, cc
zc, cc = report("C15", np.diag([a15] * 3), fr, ty, 4.0)
assert zc == [(16, 4)] and cc == [(12, 6)]
# C14 con el mismo corte relativo a la distancia Zr-Zr (a*sqrt(3)/... ~ misma escala que C15 tras escalar por a14/a15*sqrt(2)... se usa 4.35*(a14*sqrt(2)... )
zc, cc = report("C14", cell14, fr2, ty2, 4.0)
assert zc == [(16, 4)] and cc == [(12, 6)], "C14 mal construida"
print("escritos: refs/C15_init.data (648 at), refs/C14_init.data (%d at)" % len(F2))
