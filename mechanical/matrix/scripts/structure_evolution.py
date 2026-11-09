#!/usr/bin/env python3
"""Evolucion de la estructura local (Voronoi, OVITO) durante la traccion, cuadro por cuadro (SOLO LECTURA de los dumps; escribe archivos nuevos).
Por caso <caso>/*.dump (41 cuadros, eps = 0.005 k, k = 0..40) calcula con VoronoiAnalysisModifier de parametros default (sin radios ni umbrales de cara: los que reproducen la Fig. 6 de Mendez; ver structural_analysis/voronoi/voronoi_51.py):
  - CN por atomo y clase de poliedro (0 otro, 1 ico <0,0,12,0>, 2 Z14 <0,0,12,2>, 3 Z15 <0,0,12,3>, 4 Z16 <0,0,12,4>) -> <caso>/analysis/structure/voro_frames.npz  (id, type, cn[41,N] uint8, cls[41,N] uint8)
  - filas tidy (caso, muestra, T, erate, dir, eps, especie, ...) -> <caso>/analysis/structure/voro_frames.csv: media de CN, fraccion con cada CN (6..23) y fraccion de ico/Z14/Z15/Z16 por especie ('all','Zr','Cu','Ag').
Mapeo type->especie: 1 Zr, 2 Cu, 3 Ag. Los casos con voro_frames.npz mas nuevo que el dump se saltean (--force para rehacer).
Uso (modulo anaconda con OVITO): python structure_evolution.py [--sample 2e10,1e11,1e12,qgrad_ext] [--T 300] [--erate 1e8,1e9,1e10] [--only caso] [--workers 6] [--frames 0,2,...]"""
import os, sys, csv, glob, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
M = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
SP = {1: "Zr", 2: "Cu", 3: "Ag"}; IDX = {1: (0, 0, 12, 0), 2: (0, 0, 12, 2), 3: (0, 0, 12, 3), 4: (0, 0, 12, 4)}; CNR = list(range(6, 24)); NF = 41

def run_case(row):
    from ovito.io import import_file
    from ovito.modifiers import VoronoiAnalysisModifier
    p = os.path.realpath(os.path.join(M, row["path"])); dump = glob.glob(f"{p}/*.dump")[0]; out = f"{p}/analysis/structure"; os.makedirs(out, exist_ok=True)
    npz = f"{out}/voro_frames.npz"
    if "--force" not in sys.argv and os.path.exists(npz) and os.path.getmtime(npz) >= os.path.getmtime(dump): return row["case_id"], "saltado"
    pipe = import_file(dump); pipe.modifiers.append(VoronoiAnalysisModifier(compute_indices=True)); nfr = pipe.source.num_frames
    assert nfr == NF, f"{row['case_id']}: {nfr} cuadros"
    frames = [int(x) for x in arg("--frames").split(",")] if arg("--frames") else list(range(NF))
    cn = np.zeros((NF, 65536), np.uint8); cls = np.zeros((NF, 65536), np.uint8); rows = []; t0 = time.time()
    for k in frames:
        d = pipe.compute(k); ids = d.particles["Particle Identifier"][...]; o = np.argsort(ids); ids = ids[o]; t = d.particles["Particle Type"][...][o]
        c = d.particles["Coordination"][...][o]; vi = d.particles["Voronoi Index"][...][o][:, 2:6]
        cl = np.zeros(len(c), np.uint8)
        for code, ix in IDX.items(): cl[(vi == np.array(ix)).all(axis=1)] = code
        cn[k] = c; cls[k] = cl
        for lab, m in [("all", np.ones(len(t), bool))] + [(SP[a], t == a) for a in SP]:
            h = np.bincount(c[m], minlength=24)[6:24] / m.sum()
            rows.append([row["case_id"], row["sample"], row["T"], row["erate"], row["dir"], f"{0.005*k:.3f}", lab, int(m.sum()), f"{c[m].mean():.4f}"] + [f"{x:.5f}" for x in h]
                        + [f"{(cl[m] == code).mean():.5f}" for code in IDX])
    if not arg("--frames"):
        np.savez_compressed(npz, id=ids, type=t, cn=cn, cls=cls)
        with open(f"{out}/voro_frames.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["case_id", "sample", "T", "erate", "dir", "eps", "species", "N", "cn_mean"] + [f"cn{n}" for n in CNR] + ["ico", "Z14", "Z15", "Z16"]); w.writerows(rows)
    return row["case_id"], f"ok ({time.time()-t0:.0f} s, {len(frames)} cuadros)"

if __name__ == "__main__":
    sm = arg("--sample", "2e10,1e11,1e12,qgrad_ext").split(","); Ts = [float(x) for x in arg("--T", "300").split(",")]; ers = arg("--erate", "1e8,1e9,1e10").split(","); only = arg("--only")
    rows = [r for r in csv.DictReader(open(f"{M}/manifest.csv")) if r["stage"] != "extra" and r["status"].startswith("ok") and r["sample"] in sm and float(r["T"]) in Ts and r["erate"] in ers and (not only or r["case_id"] == only)]
    print(len(rows), "casos", flush=True); nw = int(arg("--workers", 6))
    if nw > 1 and len(rows) > 1:
        from multiprocessing import Pool
        with Pool(nw) as pool:
            for cid, st in pool.imap_unordered(run_case, rows): print(cid, st, flush=True)
    else:
        for r in rows: print(*run_case(r), flush=True)
