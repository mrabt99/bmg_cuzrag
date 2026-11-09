#!/usr/bin/env python3
"""Forma del mayor cluster de atomos activos (D2min > umbral relativo; vecinos <= 3.75 A, PBC) en las ventanas eps 0.12-0.20: tamano, extension por eje como fraccion de la caja
(extension = L menos el mayor hueco vacio a lo largo del eje, con PBC) y si toca/atraviesa la caja. Banda de corte: extension ~ 1 en 2 ejes y pequena en el tercero; red 3D: ~1 en los tres.
Uso (OVITO + scipy): python stz_largest_cluster_shape.py [--sample 2e10,qgrad_ext] [--erate 1e8,1e9]. Salida: matrix/stz_mayor_cluster.csv"""
import os, sys, csv, glob, warnings; warnings.filterwarnings("ignore")
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from ovito.io import import_file
M = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
def span(x, L):
    s = np.sort(x); g = np.diff(np.concatenate([s, [s[0] + L]])); return (L - g.max()) / L
rows = []
for s in arg("--sample", "2e10,1e11,1e12,qgrad_ext").split(","):
    for er in arg("--erate", "1e8,1e9,1e10").split(","):
        for d in "xyz":
            c = f"{s}_300K_{er}_{d}"; z = np.load(f"{M}/{c}/analysis/stz/stz_d2min.npz"); pipe = import_file(glob.glob(f"{M}/{c}/*.dump")[0]); ax = "xyz".index(d); out = []
            for w in range(11, 20):
                dd = pipe.compute(2 * (w + 1)); o = np.argsort(dd.particles["Particle Identifier"][...]); cell = np.array(dd.cell[...]); L = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
                pos = (dd.particles["Position"][...][o].astype(np.float64) - cell[:, 3]) % L; sel = np.flatnonzero(z["d2min"][w] > float(z["thr_rel"]))
                pr = cKDTree(pos[sel], boxsize=L).query_pairs(3.75, output_type="ndarray"); lab = connected_components(coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(sel), len(sel))), directed=False)[1]
                big = np.argmax(np.bincount(lab)); P = pos[sel][lab == big]; sp = [span(P[:, i], L[i]) for i in range(3)]
                out.append([len(P), sp[ax], *[sp[i] for i in range(3) if i != ax]])
            a = np.array(out).mean(0); rows.append([s, er, d, f"{a[0]:.0f}", f"{a[1]:.2f}", f"{a[2]:.2f}", f"{a[3]:.2f}"])
    print(s, flush=True)
with open(f"{M}/stz_mayor_cluster.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(["sample", "erate", "dir", "tam_mayor_cluster_medio", "extension_eje_traccion", "extension_lateral_1", "extension_lateral_2"]); w.writerows(rows)
