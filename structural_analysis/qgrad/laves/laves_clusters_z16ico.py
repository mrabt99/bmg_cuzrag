#!/usr/bin/env python3
"""Agrupacion espacial de los poliedros que forman una Laves: Zr con Voronoi <0,0,12,4> (Z16) y Cu con <0,0,12,0> (ico), SIN exigir
BOO ni tetraedro (criterio laxo, es la definicion 'visual' de Ryltsev: distribucion espacial de centros ico/Z16).
Clusters = componentes conexas del grafo de vecinos Voronoi restringido a esos atomos. Reporta: N atomos, tamano del cluster mayor,
su composicion (Zr/Cu/Ag) y la fraccion de Ag de sus VECINOS (donde queda el Ag). Comparacion qgrad vs 2e10 (control) y 1e11/1e12.
Salida: laves_clusters_z16ico.csv, fig_laves_cluster_mayor.png"""
import csv, os
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Z16 = np.array([0, 0, 12, 4]); ICO = np.array([0, 0, 12, 0])
T_OF = lambda lab: 300 if lab.startswith("17") else 1900 - 100 * (int(lab) - 1)
rows = []
for st, labs in {"2e10": ["05", "08", "10", "12", "13", "14", "15", "17"], "1e11": ["13", "17"], "1e12": ["13", "17"], "qgrad": ["08", "10", "12", "13", "14", "15", "16", "17r" if os.path.exists("features/qgrad/17r.npz") else "17u"]}.items():
    for lab in labs:
        F = dict(np.load(f"features/{st}/{lab}.npz")); t = F["type"]; n = len(t)
        mask = ((t == 1) & (F["vindex"] == Z16).all(1)) | ((t == 2) & (F["vindex"] == ICO).all(1))
        ptr, idx = F["neigh_ptr"], F["neigh_idx"]; src = np.repeat(np.arange(n), np.diff(ptr)); m = mask[src] & mask[idx]
        G = csr_matrix((np.ones(m.sum()), (src[m], idx[m])), shape=(n, n)); _, cl = connected_components(G, directed=False)
        sz = np.bincount(cl[mask]); big = np.argmax(sz); inb = mask & (cl == big); k = int(inb.sum())
        comp = [(t[inb] == s).sum() / k for s in (1, 2, 3)]
        # vecinos del cluster mayor (fuera de el): fraccion de Ag
        nb = np.unique(idx[np.isin(src, np.where(inb)[0])]); nb = nb[~inb[nb]]
        rows.append([st, lab, T_OF(lab), int(mask.sum()), int((sz[sz > 0] >= 2).sum()), k, k / mask.sum(), *[round(c, 3) for c in comp], round((t[nb] == 3).mean(), 3)])
with open("laves_clusters_z16ico.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["set", "config", "T", "N_Z16Zr+icoCu", "n_clusters(>=2)", "cluster_mayor", "frac_en_mayor", "x_Zr", "x_Cu", "x_Ag", "x_Ag_en_vecinos"]); w.writerows(rows)
print(f"{'set':6s} {'cfg':4s} {'T':>5s} {'N':>6s} {'clus>=2':>8s} {'mayor':>7s} {'frac':>6s}  Zr/Cu/Ag(mayor)   Ag(vecinos)")
for r in rows: print(f"{r[0]:6s} {r[1]:4s} {r[2]:5d} {r[3]:6d} {r[4]:8d} {r[5]:7d} {r[6]:6.3f}  {r[7]:.2f}/{r[8]:.2f}/{r[9]:.2f}   {r[10]:.3f}")
COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad": "#7b3fbf"}
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for st in COL:
    s = sorted([r for r in rows if r[0] == st], key=lambda r: r[2]); ax.plot([r[2] for r in s], [r[5] for r in s], "-o", ms=4, color=COL[st], label=st)
ax.set(xlabel="T (K)", ylabel="cluster mayor (atomos)", title="Cluster mayor de Zr-Z16 + Cu-ico (Voronoi, sin BOO)"); ax.grid(alpha=.25); ax.legend()
fig.tight_layout(); fig.savefig("fig_laves_cluster_mayor.png", dpi=150)
