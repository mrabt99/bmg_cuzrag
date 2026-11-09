#!/usr/bin/env python3
"""Analisis de Laves (Cu2Zr) sobre features/<set>/<label>.npz (laves_features.py). numpy+scipy+matplotlib, sin ovito.
Criterios (NINGUNO se ajusta con qgrad: umbrales fijados con las referencias y verificados en los controles 2e10/1e11/1e12):
  D  (estricto, cristal)   Zr cuya subred (Cu/Ag borrados) es diamante cubico (C15) o hexagonal (C14) segun OVITO IdentifyDiamond
  T  (local, Zr)           Zr con Voronoi <0,0,12,4> (Z16), exactamente 4 vecinos Zr y esos 4 en arreglo tetraedrico:
                           max|angulo-109.47| <= TDEV_MAX (percentil 99 de las referencias C15+C14 a 300 K, redondeado hacia arriba)
  B  (cercania en BOO)     distancia de Mahalanobis (q4,q6,W6hat) a la nube de la referencia (C15+C14 a 300 K) <= chi2(3 gl, 99 %),
                           Zr con nnn=16, Cu con nnn=12
  ZrL = T & B ;  CuL = Cu ico Voronoi <0,0,12,0> con exactamente 6 vecinos Zr & B & >=3 de esos Zr en ZrL
  Laves-like = ZrL o CuL ; clusters = componentes conexas del grafo de vecinos Voronoi restringido a atomos Laves-like.
Salidas: laves_referencias.csv, laves_conteos.csv, laves_clusters.csv, fig_laves_conteos_4curvas.png, fig_laves_clusters.png"""
import csv, glob, os
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import chi2
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

Z16 = np.array([0, 0, 12, 4]); ICO = np.array([0, 0, 12, 0])
def load(path): return dict(np.load(path))
T_OF = lambda lab: 300 if lab.startswith("17") else 1900 - 100 * (int(lab) - 1)

# ---------------- referencias
refs = {k: [load(p) for p in sorted(glob.glob(f"features/refs/{k}_300K_f*.npz"))] for k in ("C15", "C14")}
zr_feat = lambda F, sel: np.c_[F["q4_16"][sel], F["q6_16"][sel], F["w6_16"][sel]]
cu_feat = lambda F, sel: np.c_[F["q4_12"][sel], F["q6_12"][sel], F["w6_12"][sel]]
ref_rows, stats = [], {}
for sp, sel_t, fe in (("Zr", 1, zr_feat), ("Cu", 2, cu_feat)):
    X = np.vstack([fe(F, F["type"] == sel_t) for k in refs for F in refs[k]])
    mu, C = X.mean(0), np.cov(X.T) + 1e-12 * np.eye(3); stats[sp] = (mu, np.linalg.inv(C))
    for k in refs:
        Y = np.vstack([fe(F, F["type"] == sel_t) for F in refs[k]])
        ref_rows.append((k, sp, "q4,q6,W6hat", f"nnn={16 if sp=='Zr' else 12}", *[f"{v:.4f}" for v in Y.mean(0)], *[f"{v:.4f}" for v in Y.std(0)]))
tdev = np.concatenate([F["tdev"][np.isfinite(F["tdev"])] for k in refs for F in refs[k]])
TDEV_MAX = float(np.ceil(np.percentile(tdev, 99)))
print(f"referencias: tdev p50/p95/p99/max = {np.percentile(tdev,50):.1f}/{np.percentile(tdev,95):.1f}/{np.percentile(tdev,99):.1f}/{tdev.max():.1f} deg  -> TDEV_MAX = {TDEV_MAX:.0f}")
CHI = chi2.ppf(0.99, 3)
def mahal(X, sp): mu, Ci = stats[sp]; d = X - mu; return np.einsum("ij,jk,ik->i", d, Ci, d)

def classify(F):
    t = F["type"]; n = len(t); zr, cu = t == 1, t == 2
    out = {}
    out["D"] = zr & (F["diamond"] > 0)
    out["T"] = zr & np.isfinite(F["tdev"]) & (F["tdev"] <= TDEV_MAX)
    B = np.zeros(n, bool); B[zr] = mahal(zr_feat(F, zr), "Zr") <= CHI; B[cu] = mahal(cu_feat(F, cu), "Cu") <= CHI
    out["B"] = B; out["Zr_B"] = zr & B
    ZrL = out["T"] & B; out["ZrL"] = ZrL
    ico = (F["vindex"] == ICO).all(1) & cu & (F["nZr"] == 6) & B
    ptr, idx = F["neigh_ptr"], F["neigh_idx"]
    nL = np.zeros(n, int)
    if ZrL.any():
        src = np.repeat(np.arange(n), np.diff(ptr)); m = ZrL[idx]; nL = np.bincount(src[m], minlength=n)
    out["CuL"] = ico & (nL >= 3)
    out["Cu_ico6Zr_B"] = ico
    out["Laves"] = out["ZrL"] | out["CuL"]
    out["Z16"] = zr & (F["vindex"] == Z16).all(1)
    out["Z16_4Zr"] = zr & np.isfinite(F["tdev"])
    return out

def clusters(F, mask):
    n = len(mask); ptr, idx = F["neigh_ptr"], F["neigh_idx"]
    src = np.repeat(np.arange(n), np.diff(ptr)); m = mask[src] & mask[idx]
    G = csr_matrix((np.ones(m.sum()), (src[m], idx[m])), shape=(n, n))
    nc, lab = connected_components(G, directed=False)
    sizes = np.bincount(lab[mask]); sizes = sizes[sizes > 0]
    return lab, sizes

# recall en referencias (deben ser ~1 por construccion; si no, el criterio esta mal)
with open("laves_referencias.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["estructura", "especie", "descriptor", "nnn", "mean_q4", "mean_q6", "mean_w6", "std_q4", "std_q6", "std_w6"]); w.writerows(ref_rows)
for k in refs:
    r = {c: np.mean([classify(F)[c].sum() / max(1, (F["type"] == (1 if c in ("D","T","Zr_B","ZrL","Z16","Z16_4Zr") else 2)).sum()) for F in refs[k]]) for c in ("D", "T", "ZrL", "CuL")}
    print(f"recall en {k} 300K (fraccion de Zr/Cu de la referencia):", {c: round(v, 3) for c, v in r.items()})

# ---------------- conteos: qgrad + controles
rows, crow, allc = [], [], {}
sets = {"2e10": [f"{i:02d}" for i in range(1, 18)], "1e11": [f"{i:02d}" for i in range(1, 18)], "1e12": [f"{i:02d}" for i in range(1, 18)],
        "qgrad": [f"{i:02d}" for i in range(6, 17)] + ["17u"] + (["17r"] if os.path.exists("features/qgrad/17r.npz") else [])}
for st, labs in sets.items():
    for lab in labs:
        F = load(f"features/{st}/{lab}.npz"); c = classify(F); T = T_OF(lab)
        rows.append([st, lab, T] + [int(c[k].sum()) for k in ("Z16", "Z16_4Zr", "D", "T", "Zr_B", "ZrL", "CuL", "Laves")])
        if c["Laves"].any():
            _, sz = clusters(F, c["Laves"]); crow.append([st, lab, T, len(sz), int(sz.max()), int((sz >= 5).sum()), int((sz >= 10).sum()), int((sz >= 20).sum())])
        else: crow.append([st, lab, T, 0, 0, 0, 0, 0])
        allc[(st, lab)] = (F, c)
hdr = ["set", "config", "T", "N_Zr_Z16", "N_Zr_Z16_4Zr", "N_Zr_diamante(D)", "N_Zr_tetra(T)", "N_Zr_BOO(B)", "N_ZrL(T&B)", "N_CuL", "N_Laves_like"]
with open("laves_conteos.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
with open("laves_clusters.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["set", "config", "T", "n_clusters", "cluster_max", "n_>=5", "n_>=10", "n_>=20"]); w.writerows(crow)
print("\n" + " ".join(f"{h[:14]:>14s}" for h in hdr))
for r in rows:
    if r[0] in ("qgrad", "2e10") and r[2] in (1500, 1000, 800, 700, 600, 500, 300): print(" ".join(f"{str(v):>14s}" for v in r))
print("\nmax N_diamante en TODAS las configs y sets:", max(r[5] for r in rows), " | max cluster Laves-like:", max(r[4] for r in crow))
print("N Zr con tetra(T) y BOO(B):  qgrad (17u) =", [r[9] for r in rows if r[0]=="qgrad" and r[1]=="17u"], " 2e10 (300K) =", [r[9] for r in rows if r[0]=="2e10" and r[1]=="17"])

# ---------------- figura
COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad": "#7b3fbf"}
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
for j, (col, tit) in enumerate([(4, "Zr Z16 con 4 vecinos Zr"), (6, "Zr Z16 tetraedrico (T)"), (8, "Zr Laves-like (T & BOO)")]):
    for st in COL:
        s = sorted([r for r in rows if r[0] == st and not (st == "qgrad" and r[2] > 1500)], key=lambda r: r[2])
        ax[j].plot([r[2] for r in s], [r[col] for r in s], "-o", ms=3, color=COL[st], label=st)
    ax[j].set(title=tit, xlabel="T (K)", ylabel="N de atomos"); ax[j].set_yscale("symlog", linthresh=1); ax[j].grid(alpha=.25); ax[j].legend(fontsize=8)
fig.suptitle("Busqueda de Laves: conteos por T (diamante en subred Zr: 0 en todas las configuraciones)"); fig.tight_layout(); fig.savefig("fig_laves_conteos_4curvas.png", dpi=150)
