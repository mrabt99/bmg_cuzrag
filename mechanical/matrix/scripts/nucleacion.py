#!/usr/bin/env python3
"""Mapeo local de nucleacion de STZ: estructura Voronoi en eps = 0 (previa a la carga) vs deformacion local posterior, atomo por atomo.
Por caso (requiere analysis/stz/stz_d2min.npz; usa dump cuadro 0 y el cuadro del pico):
  resultado por atomo: (a) 'activo antes del pico' = D2min > umbral relativo (percentil 99.9 de la 1a ventana) en alguna ventana de d_eps = 0.01 hasta la que contiene el pico; (b) eta de von Mises ACUMULADO
  (AtomicStrain, referencia eps = 0, cutoff 3.75 A) en el cuadro del pico.
  descriptores en eps = 0: especie, CN, indice de Voronoi completo <n3,n4,n5,n6>, ico <0,0,12,0>, volumen atomico de Voronoi relativo a su especie, y versiones granulares (promedio en r < 6 A, sin el propio atomo):
  densidad de centros ico y volumen medio relativo de los vecinos.
  salidas: nucleacion_casos.csv (una fila por caso: fraccion activa, AUC de modelos logisticos con descriptores crecientes con division ESPACIAL entrenar x < L/2 / probar x >= L/2, coeficientes estandarizados del
  modelo completo, Spearman de cada descriptor vs eta_pico y vs n.o de ventanas activas, correlacion espacial en voxeles ~ 10 A entre densidad ico inicial y actividad) y nucleacion_indices.csv (riesgo relativo
  P(activo | indice, especie) / P(activo | especie) para los indices con >= 150 atomos).
Uso (anaconda: OVITO + scipy + sklearn): python nucleacion.py [--sample ...] [--T 300] [--erate ...] [--workers 6]"""
import os, sys, csv, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
M = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
RCG = 6.0; CUT = 3.75; NMIN = 150

def run_case(row):
    from ovito.io import import_file
    from ovito.modifiers import VoronoiAnalysisModifier, AtomicStrainModifier
    from ovito.pipeline import ReferenceConfigurationModifier as R
    from scipy.spatial import cKDTree
    from scipy.stats import spearmanr
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    c = row["case_id"]; p = os.path.realpath(os.path.join(M, row["path"])); dump = glob.glob(f"{p}/*.dump")[0]
    rs = [r for r in csv.DictReader(open(f"{M}/resumen.csv")) if r["case_id"] == c][0]; pk = int(np.clip(round(float(rs["eps_sigma_max"]) / 0.01) - 1, 0, 19))
    z = np.load(f"{p}/analysis/stz/stz_d2min.npz"); d2 = z["d2min"]; thr = float(z["thr_rel"]); nact = (d2[: pk + 1] > thr).sum(0); act = nact > 0
    pv = import_file(dump); pv.modifiers.append(VoronoiAnalysisModifier(compute_indices=True)); d = pv.compute(0); o = np.argsort(d.particles["Particle Identifier"][...])
    t = d.particles["Particle Type"][...][o]; cn = d.particles["Coordination"][...][o].astype(float); vol = d.particles["Atomic Volume"][...][o]; vi = d.particles["Voronoi Index"][...][o][:, 2:6]
    cell = np.array(d.cell[...]); L = np.array([cell[0, 0], cell[1, 1], cell[2, 2]]); pos = (d.particles["Position"][...][o].astype(np.float64) - cell[:, 3]) % L
    ico = (vi == np.array([0, 0, 12, 0])).all(1); idx = np.array([",".join(map(str, r)) for r in vi])
    vrel = vol.copy()
    for s in (1, 2, 3): vrel[t == s] /= vol[t == s].mean()
    tr = cKDTree(pos, boxsize=L); nb = tr.query_ball_point(pos, RCG); cg_ico = np.array([ico[[j for j in n if j != i]].mean() if len(n) > 1 else 0 for i, n in enumerate(nb)])
    cg_vol = np.array([vrel[[j for j in n if j != i]].mean() if len(n) > 1 else 1 for i, n in enumerate(nb)])
    pe = import_file(dump); pe.modifiers.append(AtomicStrainModifier(cutoff=CUT, reference_frame=0, affine_mapping=R.AffineMapping.ToReference, eliminate_cell_deformation=False))
    dd = pe.compute(2 * (pk + 1)); eta = dd.particles["Shear Strain"][...][np.argsort(dd.particles["Particle Identifier"][...])]
    X = {"especie": np.c_[t == 2, t == 3].astype(float), "CN": cn[:, None], "ico": ico[:, None].astype(float), "vol_rel": vrel[:, None], "cg_ico": cg_ico[:, None], "cg_vol": cg_vol[:, None]}
    sets = [("M0_especie", ["especie"]), ("M1_+CN", ["especie", "CN"]), ("M2_+ico", ["especie", "CN", "ico"]), ("M3_+vol", ["especie", "CN", "ico", "vol_rel"]), ("M4_+entorno", ["especie", "CN", "ico", "vol_rel", "cg_ico", "cg_vol"])]
    tr_m = pos[:, 0] < L[0] / 2; auc = {}; coef = {}
    for nm, ks in sets:
        A = np.hstack([X[k] for k in ks]); mu, sd = A.mean(0), A.std(0) + 1e-12; A = (A - mu) / sd
        m = LogisticRegression(max_iter=500).fit(A[tr_m], act[tr_m]); auc[nm] = roc_auc_score(act[~tr_m], m.predict_proba(A[~tr_m])[:, 1])
        if nm.startswith("M4"):
            names = ["Cu", "Ag", "CN", "ico", "vol_rel", "cg_ico", "cg_vol"]; coef = dict(zip(names, m.coef_[0]))
    sp = {}
    for k, v in (("CN", cn), ("ico", ico.astype(float)), ("vol_rel", vrel), ("cg_ico", cg_ico), ("cg_vol", cg_vol)):
        sp[f"rho_eta_{k}"] = spearmanr(v, eta)[0]; sp[f"rho_nact_{k}"] = spearmanr(v, nact)[0]
    nv = 10; vx = np.minimum((pos / L * nv).astype(int), nv - 1); vid = vx[:, 0] * nv * nv + vx[:, 1] * nv + vx[:, 2]; cnt = np.bincount(vid, minlength=nv ** 3)
    ok = cnt > 20; fi = np.bincount(vid, weights=ico, minlength=nv ** 3)[ok] / cnt[ok]; fa = np.bincount(vid, weights=act.astype(float), minlength=nv ** 3)[ok] / cnt[ok]; fe = np.bincount(vid, weights=eta, minlength=nv ** 3)[ok] / cnt[ok]
    vox = dict(rho_vox_ico_act=spearmanr(fi, fa)[0], rho_vox_ico_eta=spearmanr(fi, fe)[0])
    row_out = dict(case_id=c, sample=row["sample"], T=row["T"], erate=row["erate"], dir=row["dir"], eps_pico=0.01 * (pk + 1), frac_activo=act.mean(), thr=thr, **{f"AUC_{k}": v for k, v in auc.items()}, **{f"coef_{k}": v for k, v in coef.items()}, **sp, **vox)
    idx_rows = []
    for s, nm in ((1, "Zr"), (2, "Cu"), (3, "Ag")):
        ms = t == s; base = act[ms].mean()
        for u in np.unique(idx[ms]):
            m = ms & (idx == u)
            if m.sum() >= NMIN: idx_rows.append([c, row["sample"], row["T"], row["erate"], row["dir"], nm, u, int(m.sum()), f"{act[m].mean():.4f}", f"{act[m].mean()/base:.3f}", f"{eta[m].mean():.4f}"])
    return row_out, idx_rows

if __name__ == "__main__":
    sm = arg("--sample", "2e10,1e11,1e12,qgrad_ext").split(","); Ts = [float(x) for x in arg("--T", "300").split(",")]; ers = arg("--erate", "1e8,1e9,1e10").split(",")
    rows = [r for r in csv.DictReader(open(f"{M}/manifest.csv")) if r["stage"] != "extra" and r["status"].startswith("ok") and r["sample"] in sm and float(r["T"]) in Ts and r["erate"] in ers
            and os.path.exists(f"{M}/{r['path']}/analysis/stz/stz_d2min.npz")]
    print(len(rows), "casos", flush=True)
    from multiprocessing import Pool
    with Pool(int(arg("--workers", 6))) as pool: res = pool.map(run_case, rows)
    tag = arg("--tag", "")
    with open(f"{M}/nucleacion_casos{tag}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0][0].keys())); w.writeheader(); w.writerows([{k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r[0].items()} for r in res])
    with open(f"{M}/nucleacion_indices{tag}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["case_id", "sample", "T", "erate", "dir", "especie", "indice", "N", "P_activo", "riesgo_relativo", "eta_medio_pico"])
        for r in res: w.writerows(r[1])
    print("ok", flush=True)
