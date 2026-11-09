#!/usr/bin/env python3
"""Zonas de transformacion por corte (STZ) a partir del desplazamiento no afin D2min de Falk-Langer (OVITO AtomicStrainModifier, 'Nonaffine Squared Displacement', cutoff 3.75 A = primer minimo de g(r)).
Ventanas NO solapadas de d_eps = 0.01 (2 cuadros del dump, referencia = inicio de la ventana; a 1e9 s^-1 son 10 ps). Un atomo es 'activo' si D2min > D*^2.
Umbral D*^2 (relativo al ruido termico de CADA caso): percentil --pct (default 99.9) de D2min en la primera ventana (eps 0-0.01, regimen elastico: lo que ve ahi es vibracion termica). Tambien se reporta con umbrales absolutos fijos (--abs 10,20 A^2).
Agrupamiento: componentes conexas del grafo de vecinos (distancia <= 3.75 A, condiciones periodicas) entre atomos activos. CONTROL: la misma cantidad de atomos activos elegida al azar (nperm permutaciones) -> si los clusters reales
son mucho mas grandes que los del control, la actividad esta correlacionada espacialmente (STZ / percolacion), no es ruido atomo a atomo.
Salidas por caso: <caso>/analysis/stz/stz_d2min.npz (id, d2min[20,N] float32, eps_fin[20], thr) y una fila por (caso, ventana, umbral) en matrix/stz_windows.csv.
Uso (OVITO + scipy): python stz_analysis.py [--sample 2e10,1e11,1e12,qgrad_ext] [--erate 1e9] [--T 300] [--workers 6] [--pct 99.9] [--only caso]"""
import os, sys, csv, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
M = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
CUT = 3.75; NW = 20; PCT = float(arg("--pct", 99.9)); ABS = [float(x) for x in arg("--abs", "10,20").split(",")]; NPERM = int(arg("--nperm", 5))

def clusters(pos, box, sel, cut):
    from scipy.spatial import cKDTree
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    ii = np.flatnonzero(sel)
    if len(ii) == 0: return np.array([], int)
    t = cKDTree(pos[ii], boxsize=box); pr = t.query_pairs(cut, output_type="ndarray")
    g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(ii), len(ii)))
    _, lab = connected_components(g, directed=False); return np.bincount(lab)

def stats(sz, n):
    if len(sz) == 0: return dict(n_cl=0, n_cl5=0, largest=0, mean_sz_ge2=0.0, frac_in_ge5=0.0)
    return dict(n_cl=int(len(sz)), n_cl5=int((sz >= 5).sum()), largest=int(sz.max()), mean_sz_ge2=float(sz[sz >= 2].mean()) if (sz >= 2).any() else 0.0, frac_in_ge5=float(sz[sz >= 5].sum() / n) if n else 0.0)

def run_case(row):
    from ovito.io import import_file
    from ovito.modifiers import AtomicStrainModifier
    p = os.path.realpath(os.path.join(M, row["path"])); dump = glob.glob(f"{p}/*.dump")[0]; out = f"{p}/analysis/stz"; os.makedirs(out, exist_ok=True)
    pipe = import_file(dump); pipe.modifiers.append(AtomicStrainModifier(cutoff=CUT, output_nonaffine_squared_displacements=True, use_frame_offset=True, frame_offset=-2, affine_mapping=AtomicStrainModifier.AffineMapping.ToReference))
    D2 = np.zeros((NW, 65536), np.float32); POS = []; BOX = []; EF = []
    for w in range(NW):
        k = 2 * (w + 1); d = pipe.compute(k); ids = d.particles["Particle Identifier"][...]; o = np.argsort(ids)
        D2[w] = d.particles["Nonaffine Squared Displacement"][...][o]; c = np.array(d.cell[...]); L = np.array([c[0, 0], c[1, 1], c[2, 2]]); org = c[:, 3]
        POS.append((d.particles["Position"][...][o].astype(np.float64) - org) % L); BOX.append(L); EF.append(0.005 * k)
    thr_rel = float(np.percentile(D2[0], PCT)); rows = []; rng = np.random.default_rng(1)
    np.savez_compressed(f"{out}/stz_d2min.npz", id=np.sort(ids), d2min=D2, eps_fin=np.array(EF), thr_rel=thr_rel)
    for lab, thr in [(f"p{PCT}_1a_ventana", thr_rel)] + [(f"abs{a:g}", a) for a in ABS]:
        for w in range(NW):
            sel = D2[w] > thr; n = int(sel.sum()); s = stats(clusters(POS[w], BOX[w], sel, CUT), n)
            nu = [stats(clusters(POS[w], BOX[w], rng.permutation(sel), CUT), n) for _ in range(NPERM)]
            nm = {k: float(np.mean([x[k] for x in nu])) for k in s}
            rows.append([row["case_id"], row["sample"], row["erate"], row["dir"], lab, f"{thr:.3f}", f"{EF[w]-0.01:.3f}", f"{EF[w]:.3f}", n, f"{n/65536:.5f}", f"{D2[w].mean():.4f}"]
                        + [s[k] for k in ("n_cl", "n_cl5", "largest", "mean_sz_ge2", "frac_in_ge5")] + [f"{nm[k]:.4g}" for k in ("n_cl", "n_cl5", "largest", "mean_sz_ge2", "frac_in_ge5")])
    return rows

HDR = ["case_id", "sample", "erate", "dir", "umbral", "thr_A2", "eps_ini", "eps_fin", "n_activos", "frac_activos", "D2min_medio", "n_cl", "n_cl_ge5", "mayor", "tam_medio_ge2", "frac_en_cl_ge5", "null_n_cl", "null_n_cl_ge5", "null_mayor", "null_tam_medio_ge2", "null_frac_en_cl_ge5"]
if __name__ == "__main__":
    sm = arg("--sample", "2e10,1e11,1e12,qgrad_ext").split(","); ers = arg("--erate", "1e9").split(","); Ts = [float(x) for x in arg("--T", "300").split(",")]; only = arg("--only")
    rows = [r for r in csv.DictReader(open(f"{M}/manifest.csv")) if r["stage"] != "extra" and r["status"].startswith("ok") and r["sample"] in sm and r["erate"] in ers and float(r["T"]) in Ts and (not only or r["case_id"] == only)]
    print(len(rows), "casos", flush=True)
    from multiprocessing import Pool
    outp = f"{M}/stz_windows.csv"; prev = []
    if os.path.exists(outp): prev = [r for r in csv.reader(open(outp))][1:]; done = {r["case_id"] for r in rows}; prev = [r for r in prev if r[0] not in done]
    with Pool(int(arg("--workers", 6))) as pool: res = pool.map(run_case, rows)
    with open(outp, "w", newline="") as f:
        w = csv.writer(f); w.writerow(HDR); w.writerows(prev)
        for r in res: w.writerows(r)
    print("->", outp)
