#!/usr/bin/env python3
"""Limite de deteccion: implanta esferas de radio R de cristal Cu2Zr (C15 o C14, termalizado a 300 K, refs/*_300K.pos.dump ultimo frame)
en el vidrio 2e10 a 300 K (snaps/17.dump), reemplazando los atomos de vidrio a menos de R+1.2 A del centro (sin relajar).
Mide con los detectores de laves_features (D = diamante en subred Zr; T = Zr-Z16 con 4 Zr tetraedricos) y con S(k).
Corre con anaconda (OVITO). Salida: implant/<tipo>_R<R>.dump, laves_implant_deteccion.csv."""
import os, sys, csv, warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import numpy as np
import laves_features as LF
SA = "/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis/"
def read_dump(path, last=True):
    frames, i = [], 0
    L = open(path).read().split("ITEM: TIMESTEP")[1:]
    txt = L[-1] if last else L[0]; lines = txt.strip().split("\n")
    n = int(lines[2]); hdr = lines[3]; b = [list(map(float, l.split())) for l in lines[4:7]]
    atoms = np.array([l.split() for l in lines[8:8 + n]], float); atoms = atoms[np.argsort(atoms[:, 0])]
    return atoms, np.array(b), hdr
glass, gb, _ = read_dump(SA + "2e10/snaps/17.dump"); gL = gb[:, 1] - gb[:, 0]; assert np.allclose(gL, gL[0])
center = gb[:, 0] + gL / 2
def crystal_fragment(tag):
    at, b, hdr = read_dump(f"refs/{tag}_300K.pos.dump")
    if "xy" in hdr:   # triclinico: b = [xlo_bound xhi_bound xy]...
        xy, xz, yz = b[0, 2], b[1, 2], b[2, 2]
        xlo = b[0, 0] - min(0, xy, xz, xy + xz); xhi = b[0, 1] - max(0, xy, xz, xy + xz); ylo = b[1, 0] - min(0, yz); yhi = b[1, 1] - max(0, yz)
        A = np.array([[xhi - xlo, xy, xz], [0, yhi - ylo, yz], [0, 0, b[2, 1] - b[2, 0]]]).T; o = np.array([xlo, ylo, b[2, 0]])
    else:
        A = np.diag(b[:, 1] - b[:, 0]); o = b[:, 0]
    reps = range(-3, 4); P, T = [], []
    for i in reps:
        for j in reps:
            for k in reps:
                P.append(at[:, 2:5] + (A @ np.array([i, j, k]))[None, :]); T.append(at[:, 1])
    return np.vstack(P), np.concatenate(T), A
def implant(tag, R):
    P, T, A = crystal_fragment(tag); c0 = P[np.argmin(np.linalg.norm(P - P.mean(0), axis=1))]
    sel = np.linalg.norm(P - c0, axis=1) <= R
    Pc, Tc = P[sel] - c0 + center, T[sel]
    dg = np.linalg.norm(glass[:, 2:5] - center, axis=1); keep = dg > R + 1.2
    atoms = np.vstack([np.c_[glass[keep, 1], glass[keep, 2:5]], np.c_[Tc, Pc]])
    out = f"implant/{tag}_R{R:g}.dump"
    with open(out, "w") as f:
        f.write(f"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{len(atoms)}\nITEM: BOX BOUNDS pp pp pp\n")
        for d in range(3): f.write(f"{gb[d,0]:.8f} {gb[d,1]:.8f}\n")
        f.write("ITEM: ATOMS id type x y z\n")
        for i, a in enumerate(atoms, 1): f.write(f"{i} {int(a[0])} {a[1]:.6f} {a[2]:.6f} {a[3]:.6f}\n")
    return out, int(sel.sum()), int(keep.sum())
if __name__ == "__main__":
    os.makedirs("implant", exist_ok=True); rows = []
    for tag in ("C15", "C14"):
        for R in (8, 10, 12, 15, 20):
            out, nc, ng = implant(tag, R)
            F = LF.compute(out, None, 0); t = F["type"]
            nz = int(((t == 1)).sum()); D = int((F["diamond"] > 0).sum()); Tt = int((np.isfinite(F["tdev"]) & (F["tdev"] <= 10)).sum())
            zr_in = int((Tc := None) is None) and int(0)
            rows.append([tag, R, nc, D, Tt]); print(f"{tag} R={R:2d} A: atomos del cristal={nc:5d}  -> D (diamante Zr)={D:4d}  T (tetra Z16)={Tt:4d}", flush=True)
    with open("laves_implant_deteccion.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["cristal", "R_A", "N_atomos_cristal", "N_Zr_diamante_detectados", "N_Zr_tetra_detectados"]); w.writerows(rows)
