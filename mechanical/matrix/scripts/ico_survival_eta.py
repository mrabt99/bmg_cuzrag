#!/usr/bin/env python3
"""Replica del analisis de supervivencia de icosaedros de Mendez & Jaroszewicz (J. Non-Cryst. Solids 693, 124394, Sec. 2.4 y 3.4.3): centros ico <0,0,12,0> (Voronoi) identificados en eps = 0, seguidos atomo por atomo;
en cada cuadro, eta de von Mises de CADA atomo respecto de la MISMA configuracion de referencia (eps = 0; AtomicStrainModifier, acumulado, no incremental); atomo 'elastico' si eta < eta* = 0.05, 'activo' si eta >= eta*.
Se reporta la probabilidad de supervivencia del ico inicial en cada clase. El paper da a eps ~ 0.01: 73 % (elasticos, n = 1484) vs 57 % (activos, n = 979).
Parametros no dados en el paper y por eso barridos aca: cutoff del AtomicStrain (default de OVITO 3.0 A; 3.75 A = primer minimo de g(r)). Se agrega un CONTROL de ruido termico: el mismo corte en eta aplicado a eps = 0.005 y a un cuadro (eps = 0.01) pero
mirando tambien cuanto del 'activo' se debe a la amplitud vibracional (eta mediano de los atomos NO ico).
Usa Voronoi ya calculado (<caso>/analysis/structure/voro_frames.npz). Salida: matrix/ico_survival_eta.csv. Uso (anaconda/OVITO): python ico_survival_eta.py [--erate 1e9] [--workers 6]"""
import os, sys, csv, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
M = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
FR = [1, 2, 4, 8, 12, 20, 40]; CUT = [3.0, 3.75]; ETA_STAR = 0.05

def run_case(row):
    from ovito.io import import_file
    from ovito.modifiers import AtomicStrainModifier
    from ovito.pipeline import ReferenceConfigurationModifier as R
    p = os.path.realpath(os.path.join(M, row["path"])); dump = glob.glob(f"{p}/*.dump")[0]; z = np.load(f"{p}/analysis/structure/voro_frames.npz")
    ico0 = z["cls"][0] == 1; typ = z["type"]; out = []
    for cut in CUT:
        pipe = import_file(dump); pipe.modifiers.append(AtomicStrainModifier(cutoff=cut, reference_frame=0, affine_mapping=R.AffineMapping.ToReference, eliminate_cell_deformation=False))
        for k in FR:
            d = pipe.compute(k); ids = d.particles["Particle Identifier"][...]; o = np.argsort(ids); eta = d.particles["Shear Strain"][...][o]
            surv = z["cls"][k] == 1; act = eta >= ETA_STAR
            for lab, m in (("ico_elastico", ico0 & ~act), ("ico_activo", ico0 & act)):
                out.append([row["case_id"], row["sample"], row["dir"], cut, f"{0.005*k:.3f}", lab, int(m.sum()), f"{surv[m].mean():.4f}", f"{np.median(eta[ico0]):.4f}", f"{np.median(eta[~ico0]):.4f}", f"{act.mean():.4f}"])
    return out

if __name__ == "__main__":
    er = arg("--erate", "1e9"); sm = arg("--sample", "2e10,1e11,1e12,qgrad_ext").split(",")
    rows = [r for r in csv.DictReader(open(f"{M}/manifest.csv")) if r["stage"] != "extra" and r["status"].startswith("ok") and r["sample"] in sm and float(r["T"]) == 300 and r["erate"] == er]
    from multiprocessing import Pool
    with Pool(int(arg("--workers", 6))) as pool: res = pool.map(run_case, rows)
    with open(f"{M}/ico_survival_eta_{er}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["case_id", "sample", "dir", "cutoff", "eps", "clase", "n", "supervivencia", "eta_med_ico", "eta_med_resto", "frac_activos_total"])
        for r in res: w.writerows(r)
    print(len(rows), "casos ->", f"{M}/ico_survival_eta_{er}.csv")
