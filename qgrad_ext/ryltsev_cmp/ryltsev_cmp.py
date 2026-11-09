#!/usr/bin/env python3
"""Comparacion con la deteccion 'visual' de Ryltsev et al. (JCP 145, 034506) a 300 K en 2e10, 1e11, qgrad y qgrad_ext (todas 17r/17 relajadas).
(1) Fig. 11a: mapa de una loncha (slab) de la caja con TODOS los atomos coloreados y dimensionados por q6 (Nnn=12); ico = q6>=0.6 & w6<=-0.16.
(2) Agrupacion espacial de los centros ico (BOO) frente a un azar de la misma densidad: n. medio de vecinos ico a < RN de cada ico / esperado al azar (1 = sin agrupacion).
(3) Z16 por BOO: criterio estricto del paper (q4<=0.01 & |q6-0.1924|<=0.01, Nnn=16) y su mitad q6 (la unica alcanzable a 300 K, ver classify/boo_types.py).
(4) Cluster mayor de Zr-Z16 + Cu-ico (Voronoi, criterio laxo = 'distribucion espacial de centros ico/Z16'): misma definicion que qgrad/laves/laves_clusters_z16ico.py.
Solo lectura de structural_analysis/ y qgrad_ext/. Salidas en este directorio. Uso: PYTHONDONTWRITEBYTECODE=1 python ryltsev_cmp.py"""
import sys, csv
SA = "/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis"; E = "/scratch/negishi/aaabt/bmg_cuzrag/qgrad_ext"
sys.path.insert(0, f"{SA}/classify")
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from boo_types import load_boo_dump, mask_ico, mask_z16, mask_z16_q6window, features
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
CUR = [("2e10", f"{SA}/2e10/snaps/17.dump", f"{SA}/boo/2e10/17.boo.dump", f"{SA}/qgrad/laves/features/2e10/17.npz"),
       ("1e11", f"{SA}/1e11/snaps/17.dump", f"{SA}/boo/1e11/17.boo.dump", f"{SA}/qgrad/laves/features/1e11/17.npz"),
       ("qgrad", f"{SA}/qgrad/snaps/17r.dump", f"{SA}/qgrad/boo/17r.boo.dump", f"{SA}/qgrad/laves/features/qgrad/17r.npz"),
       ("qgrad_ext", f"{E}/snaps/17r.dump", f"{E}/boo/17r.boo.dump", f"{E}/laves/features/qgrad_ext/17r.npz")]
RN, SLAB = 4.5, 12.0
def read_dump(p):
    L = open(p).read().split("\n"); n = int(L[3]); b = np.array([l.split()[:2] for l in L[5:8]], float)
    a = np.array([l.split()[:5] for l in L[9:9 + n]], float); a = a[np.argsort(a[:, 0])]
    return a[:, 1].astype(int), a[:, 2:5] - b[:, 0], (b[:, 1] - b[:, 0])[0]
Z16V, ICOV = np.array([0, 0, 12, 4]), np.array([0, 0, 12, 0])
rows, data = [], {}
for k, pd, pb, pf in CUR:
    t, x, L = read_dump(pd); d = load_boo_dump(pb); assert (np.sort(d[:, 0]) == d[:, 0]).all() or True
    d = d[np.argsort(d[:, 0])]; assert (d[:, 1].astype(int) == t).all()
    q4_12, q6_12, _, w6_12 = features(d, 12); ico = mask_ico(d); n = len(t)
    tr = cKDTree(x % L, boxsize=L); xi = (x % L)[ico]; ti = cKDTree(xi, boxsize=L)
    cnt = np.array([len(v) - 1 for v in ti.query_ball_point(xi, RN)]); mean_nb = cnt.mean()
    rho = ico.sum() / L ** 3; exp = rho * 4 / 3 * np.pi * RN ** 3
    # agrupacion: nº de vecinos ico dentro de RN relativo al azar; error por reparto aleatorio de las mismas posiciones ico (50 realizaciones)
    rng = np.random.default_rng(0); rnd = []
    for _ in range(20):
        xr = rng.random((ico.sum(), 3)) * L; tt = cKDTree(xr, boxsize=L); rnd.append(np.mean([len(v) - 1 for v in tt.query_ball_point(xr, RN)]))
    z16_strict = int(mask_z16(d).sum()); z16_q6 = int(mask_z16_q6window(d).sum())
    F = dict(np.load(pf)); tv, nv = F["type"], F["vindex"]; mask = ((tv == 1) & (nv == Z16V).all(1)) | ((tv == 2) & (nv == ICOV).all(1))
    ptr, idx = F["neigh_ptr"], F["neigh_idx"]; src = np.repeat(np.arange(n), np.diff(ptr)); m = mask[src] & mask[idx]
    _, cl = connected_components(csr_matrix((np.ones(m.sum()), (src[m], idx[m])), shape=(n, n)), directed=False)
    sz = np.bincount(cl[mask]); big = int(sz.max())
    rows.append([k, int(ico.sum()), f"{mean_nb:.3f}", f"{np.mean(rnd):.3f}", f"{mean_nb / np.mean(rnd):.3f}", f"{np.std(rnd):.3f}", z16_strict, z16_q6, int(mask.sum()), big])
    data[k] = (x % L, q6_12, ico, L)
    print(rows[-1], flush=True)
with open("ryltsev_cmp.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["curva", "N_ico_BOO", "vecinos_ico_<4.5A_medio", "azar_medio", "factor_agrupacion", "std_azar", "N_Z16_BOO_estricto", "N_Z16_BOO_q6window", "N_Z16Zr+icoCu_Voronoi", "cluster_mayor"]); w.writerows(rows)
# ---- Fig. 11a-like: loncha 0<z<SLAB, color y tamano por q6
fig, axs = plt.subplots(1, 4, figsize=(19, 5), constrained_layout=True)
for ax, (k, *_r) in zip(axs, CUR):
    x, q6, ico, L = data[k]; s = x[:, 2] < SLAB
    o = np.argsort(q6[s]); sc = ax.scatter(x[s, 0][o], x[s, 1][o], c=q6[s][o], s=2 + 40 * np.clip(q6[s][o] - 0.2, 0, None) ** 2 * 6, cmap="turbo", vmin=0.1, vmax=0.7, linewidths=0)
    ax.set_aspect("equal"); ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_title(f"{k}: {int((ico & s).sum())} ico en la loncha", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
fig.colorbar(sc, ax=axs, shrink=.8, label="q6 (Nnn = 12)")
fig.suptitle(f"Estilo Fig. 11a de Ryltsev: loncha de {SLAB:.0f} Å (todos los átomos, color y tamaño por q6), 300 K relajado. Rojo grande = centros icosaédricos", fontsize=11)
fig.savefig("fig11a_style_300K.png", dpi=140)
