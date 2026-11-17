#!/usr/bin/env python3
"""¿Los icosaedros (Cu) y los Frank-Kasper (Zr) forman clusters? Vidrios a 300 K: 1e12 (UNIFICADA), 1e11, 2e10 y qgrad_ext.
Rehace la idea de ryltsev_cmp.py (17/11/2026) con las muestras correctas: sin qgrad (= 2e10) y con la 1e12 unificada en vez
de nada. ryltsev_cmp.py y sus salidas no se tocan.

Para cada vidrio y cada tipo de poliedro:
  * clusters = componentes conexas del grafo de VECINOS DE VORONOI (OVITO default; features/<set>/<lab>.npz de
    laves_features.py) restringido a los átomos de ese tipo;
  * estadísticas: N, fracción de esos átomos que está en un cluster de 2 o más, tamaño medio de los clusters de 2 o más, cluster mayor;
  * REFERENCIA AL AZAR: las mismas etiquetas repartidas al azar entre los átomos de la misma especie (20 repeticiones).
    Observado/azar > 1 = se agrupan más que por casualidad;
  * además, para ico BOO, el factor de agrupación espacial de ryltsev_cmp.py (vecinos ico a < 4,5 Å / azar en el espacio).
Tipos: Cu ico Voronoi <0,0,12,0>; Cu ico BOO (q6>=0,6 y W6<=-0,16, Nb=12); Zr FK = <0,0,12,2>+<0,0,12,3>+<0,0,12,4>;
       Zr Z16 <0,0,12,4>; y la mezcla de Ryltsev Zr-Z16 + Cu-ico (Voronoi).
Salida: clusters_4tasas.csv. Uso: PYTHONDONTWRITEBYTECODE=1 /scratch/negishi/aaabt/venv_structural/bin/python clusters_4tasas.py
"""
import csv
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

SA = "/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis"; E = "/scratch/negishi/aaabt/bmg_cuzrag/qgrad_ext"
MUESTRAS = [("1e12", f"{SA}/1e12u/snaps/17.dump", f"{SA}/qgrad/laves/features/1e12u/17.npz"),
            ("1e11", f"{SA}/1e11/snaps/17.dump", f"{SA}/qgrad/laves/features/1e11/17.npz"),
            ("2e10", f"{SA}/2e10/snaps/17.dump", f"{SA}/qgrad/laves/features/2e10/17.npz"),
            ("qgrad_ext", f"{E}/snaps/17r.dump", f"{E}/laves/features/qgrad_ext/17r.npz")]
RN, NREP = 4.5, 20
ICO, Z14, Z15, Z16 = (np.array(v) for v in ([0, 0, 12, 0], [0, 0, 12, 2], [0, 0, 12, 3], [0, 0, 12, 4]))


def read_dump(p):
    L = open(p).read().split("\n"); n = int(L[3]); b = np.array([l.split()[:2] for l in L[5:8]], float)
    a = np.array([l.split()[:5] for l in L[9:9 + n]], float); a = a[np.argsort(a[:, 0])]
    return a[:, 1].astype(int), (a[:, 2:5] - b[:, 0]) % (b[0, 1] - b[0, 0]), b[0, 1] - b[0, 0]


def stats(mask, src, idx, n):
    m = mask[src] & mask[idx]
    _, cl = connected_components(csr_matrix((np.ones(m.sum()), (src[m], idx[m])), shape=(n, n)), directed=False)
    sz = np.bincount(cl[mask]); sz = sz[sz > 0]; ge2 = sz[sz >= 2]
    return dict(N=int(mask.sum()), frac_en_cl=ge2.sum() / mask.sum(), tam_medio=ge2.mean() if len(ge2) else 1.0, mayor=int(sz.max()))


rng = np.random.default_rng(12345)
filas = []
for nombre, pdump, pfeat in MUESTRAS:
    t, x, L = read_dump(pdump)
    F = dict(np.load(pfeat)); assert (F["type"] == t).all(), nombre
    n = len(t); v = F["vindex"]
    src = np.repeat(np.arange(n), np.diff(F["neigh_ptr"])); idx = F["neigh_idx"]
    cu, zr = t == 2, t == 1
    tipos = {"Cu ico (Voronoi)": cu & (v == ICO).all(1),
             "Cu ico (BOO)": cu & (F["q6_12"] >= 0.6) & (F["w6_12"] <= -0.16),
             "Zr FK (Z14+Z15+Z16)": zr & ((v == Z14).all(1) | (v == Z15).all(1) | (v == Z16).all(1)),
             "Zr Z16": zr & (v == Z16).all(1),
             "Zr-Z16 + Cu-ico": (zr & (v == Z16).all(1)) | (cu & (v == ICO).all(1))}
    for tipo, mask in tipos.items():
        obs = stats(mask, src, idx, n)
        nul = []
        for _ in range(NREP):
            mr = np.zeros(n, bool)
            for esp in (1, 2, 3):
                sel = np.where(t == esp)[0]; k = int((mask & (t == esp)).sum())
                if k: mr[rng.choice(sel, k, replace=False)] = True
            nul.append(stats(mr, src, idx, n))
        r = dict(muestra=nombre, tipo=tipo, N=obs["N"])
        for c in ("frac_en_cl", "tam_medio", "mayor"):
            a = np.array([d[c] for d in nul], float)
            r[c] = obs[c]; r[c + "_azar"] = a.mean(); r[c + "_azar_sd"] = a.std(); r[c + "_razon"] = obs[c] / a.mean()
        if tipo == "Cu ico (BOO)":   # factor espacial de ryltsev_cmp.py
            xi = x[mask]; ti = cKDTree(xi, boxsize=L)
            nb = np.mean([len(q) - 1 for q in ti.query_ball_point(xi, RN)])
            rnd = []
            for _ in range(NREP):
                xr = rng.random((mask.sum(), 3)) * L; tr = cKDTree(xr, boxsize=L)
                rnd.append(np.mean([len(q) - 1 for q in tr.query_ball_point(xr, RN)]))
            r["factor_espacial"] = nb / np.mean(rnd)
        filas.append(r)
        print(f"{nombre:9s} {tipo:20s} N={obs['N']:5d} | en clusters {100*obs['frac_en_cl']:5.1f}% (azar {100*r['frac_en_cl_azar']:5.1f}) "
              f"| tam medio {obs['tam_medio']:5.2f} (azar {r['tam_medio_azar']:4.2f}) | mayor {obs['mayor']:4d} "
              f"(azar {r['mayor_azar']:5.1f}±{r['mayor_azar_sd']:4.1f})" + (f" | factor espacial {r['factor_espacial']:.2f}" if "factor_espacial" in r else ""))
cols = ["muestra", "tipo", "N"] + [f"{c}{s}" for c in ("frac_en_cl", "tam_medio", "mayor") for s in ("", "_azar", "_azar_sd", "_razon")] + ["factor_espacial"]
with open("clusters_4tasas.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in filas:
        w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
