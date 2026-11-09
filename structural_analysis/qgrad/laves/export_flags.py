#!/usr/bin/env python3
"""Exporta un dump para OVITO (local) con flags por atomo: D (Zr diamante), T (Zr-Z16 tetraedrico), ZrL, CuL, Z16 (Zr Voronoi <0,0,12,4>),
CuIco (Cu Voronoi <0,0,12,0>), BigCl (1 = cluster mayor de Zr-Z16 + Cu-ico). Ejecutar con venv_structural. Uso: python export_flags.py <set> <label> <snap_dump>"""
import sys
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
st, lab, snap = sys.argv[1:4]
F = dict(np.load(f"features/{st}/{lab}.npz")); t = F["type"]; n = len(t)
L = open(snap).read().split("\n"); nb = np.array([l.split()[:2] for l in L[5:8]], float); at = np.array([l.split() for l in L[9:9 + n]], float); at = at[np.argsort(at[:, 0])]
assert np.array_equal(at[:, 1].astype(int), t)
Z16 = (t == 1) & (F["vindex"] == [0, 0, 12, 4]).all(1); ICO = (t == 2) & (F["vindex"] == [0, 0, 12, 0]).all(1); mask = Z16 | ICO
ptr, idx = F["neigh_ptr"], F["neigh_idx"]; src = np.repeat(np.arange(n), np.diff(ptr)); m = mask[src] & mask[idx]
_, cl = connected_components(csr_matrix((np.ones(m.sum()), (src[m], idx[m])), shape=(n, n)), directed=False)
big = mask & (cl == np.argmax(np.bincount(cl[mask])))
D = (F["diamond"] > 0); T = np.isfinite(F["tdev"]) & (F["tdev"] <= 10)
out = f"candidatos/{st}_{lab}_flags.dump"
with open(out, "w") as f:
    f.write(f"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp pp\n" + "".join(f"{a:.8f} {b:.8f}\n" for a, b in nb) + "ITEM: ATOMS id type x y z D T Z16 CuIco BigCl\n")
    for i in range(n): f.write(f"{int(at[i,0])} {int(at[i,1])} {at[i,2]:.5f} {at[i,3]:.5f} {at[i,4]:.5f} {int(D[i])} {int(T[i])} {int(Z16[i])} {int(ICO[i])} {int(big[i])}\n")
print(f"{out}: D={int(D.sum())} T={int(T.sum())} Z16={int(Z16.sum())} CuIco={int(ICO.sum())} cluster_mayor={int(big.sum())}")
