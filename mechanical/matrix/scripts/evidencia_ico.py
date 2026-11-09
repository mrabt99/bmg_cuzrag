#!/usr/bin/env python3
"""Tres pruebas para reforzar (o no) la idea 'mas icosaedros -> nucleacion de STZ retrasada / vidrio mas resistente' (sin simulaciones nuevas):
 (A) Control qgrad: 5 vidrios a 300 K. Se ajusta sigma_max (y delta_sigma) contra Ep inicial y contra f_ico(Cu) inicial usando los 4 vidrios SIN qgrad y se predice qgrad (por tasa de deformacion).
 (B) Pares emparejados: dentro de los Cu con CN = 12, icosaedros <0,0,12,0> vs no-icosaedros, estratificando por volumen de Voronoi (quintiles dentro del caso) y por volumen medio de los vecinos en r < 6 A (terciles);
     diferencia ponderada de P(activo antes del pico) = P(ico) - P(no-ico), por caso -> media y desvio sobre x, y, z. 'Activo' = D2min > umbral relativo (percentil 99.9 de la 1a ventana) en alguna ventana hasta el pico.
 (C) Inicio de la actividad con barras de error: eps en que la fraccion activa por ventana (stz_windows.csv, umbral relativo) cruza 0.5 % y 1 % (interpolacion lineal), por muestra y tasa (media y desvio sobre x, y, z), comparado con eps(sigma_max).
Salidas: matrix/evidencia_ico.csv (B y C por caso), evidencia_ico_resumen.txt, figs/evidencia_ico.png. Uso (anaconda + OVITO): python evidencia_ico.py"""
import os, csv, glob, warnings
warnings.filterwarnings("ignore")
from collections import defaultdict
import numpy as np
M = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = {r["case_id"]: r for r in csv.DictReader(open(f"{M}/resumen.csv"))}
SAMPLES = ["2e10", "1e11", "1e12", "qgrad_ext", "qgrad"]; ER = ["1e8", "1e9", "1e10"]

def pares(row):
    from ovito.io import import_file
    from ovito.modifiers import VoronoiAnalysisModifier
    from scipy.spatial import cKDTree
    c = row["case_id"]; p = os.path.realpath(os.path.join(M, row["path"])); dump = glob.glob(f"{p}/*.dump")[0]
    pk = int(np.clip(round(float(RES[c]["eps_sigma_max"]) / 0.01) - 1, 0, 19)); z = np.load(f"{p}/analysis/stz/stz_d2min.npz"); act = (z["d2min"][: pk + 1] > float(z["thr_rel"])).any(0)
    pv = import_file(dump); pv.modifiers.append(VoronoiAnalysisModifier(compute_indices=True)); d = pv.compute(0); o = np.argsort(d.particles["Particle Identifier"][...])
    t = d.particles["Particle Type"][...][o]; cn = d.particles["Coordination"][...][o]; vol = d.particles["Atomic Volume"][...][o]; vi = d.particles["Voronoi Index"][...][o][:, 2:6]
    ico = (vi == np.array([0, 0, 12, 0])).all(1); cell = np.array(d.cell[...]); L = np.array([cell[0, 0], cell[1, 1], cell[2, 2]]); pos = (d.particles["Position"][...][o].astype(np.float64) - cell[:, 3]) % L
    vrel = vol / vol[t == 2].mean(); nb = cKDTree(pos, boxsize=L).query_ball_point(pos, 6.0); cgv = np.array([vrel[[j for j in n if j != i]].mean() for i, n in enumerate(nb)])
    m = (t == 2) & (cn == 12); iq = np.digitize(vrel, np.quantile(vrel[m], [.2, .4, .6, .8])); ig = np.digitize(cgv, np.quantile(cgv[m], [1 / 3, 2 / 3]))
    num = 0.0; w = 0; pi_, pn_ = [], []
    for a in range(5):
        for b in range(3):
            s = m & (iq == a) & (ig == b); i1 = s & ico; i0 = s & ~ico
            if i1.sum() >= 5 and i0.sum() >= 5: num += (act[i1].mean() - act[i0].mean()) * i1.sum(); w += i1.sum(); pi_.append(act[i1].mean() * i1.sum()); pn_.append(act[i0].mean() * i1.sum())
    return dict(case_id=c, sample=row["sample"], T=row["T"], erate=row["erate"], dir=row["dir"], n_ico_pareados=int(w), P_ico=sum(pi_) / w, P_noico_CN12=sum(pn_) / w, dif=num / w, razon=(sum(pi_) / max(sum(pn_), 1e-9)))

def onset(rows, thr):
    e = np.array([float(r["eps_fin"]) for r in rows]); f = np.array([float(r["frac_activos"]) for r in rows]); o = np.argsort(e); e, f = e[o], f[o]
    for i in range(1, len(e)):
        if f[i] >= thr and f[i - 1] < thr: return e[i - 1] + (thr - f[i - 1]) / (f[i] - f[i - 1]) * (e[i] - e[i - 1])
    return np.nan

if __name__ == "__main__":
    from multiprocessing import Pool
    rows = [r for r in csv.DictReader(open(f"{M}/manifest.csv")) if r["stage"] != "extra" and r["status"].startswith("ok") and os.path.exists(f"{M}/{r['path']}/analysis/stz/stz_d2min.npz")]
    with Pool(6) as pool: P = pool.map(pares, rows)
    W = defaultdict(list)
    for r in csv.DictReader(open(f"{M}/stz_windows.csv")):
        if r["umbral"].startswith("p"): W[r["case_id"]].append(r)
    for r in P:
        c = r["case_id"]; r["onset_0.5pct"] = onset(W[c], 0.005); r["onset_1pct"] = onset(W[c], 0.01); r["eps_pico"] = float(RES[c]["eps_sigma_max"]); r["sigma_max"] = float(RES[c]["sigma_max_GPa"]); r["delta_sigma"] = float(RES[c]["delta_sigma_GPa"])
        cv = np.loadtxt(glob.glob(f"{M}/{c}/*.curve.dat")[0], comments="#"); r["PE0"] = cv[:5, 6].mean()
        v = np.load(f"{M}/{c}/analysis/structure/voro_frames.npz"); r["fico_Cu0"] = (v["cls"][0][v["type"] == 2] == 1).mean()
    with open(f"{M}/evidencia_ico.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(P[0].keys())); w.writeheader(); w.writerows([{k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()} for r in P])
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad_ext": "#b5338a", "qgrad": "#8a8985"}
    g = defaultdict(list)
    for r in P:
        if float(r["T"]) == 300: g[(r["sample"], r["erate"])].append(r)
    mean = lambda s, er, k: np.mean([r[k] for r in g[(s, er)]]); sd = lambda s, er, k: np.std([r[k] for r in g[(s, er)]], ddof=1)
    L = ["== (A) control exploratorio qgrad (NO concluyente: un vidrio por muestra y recta con 4 puntos): ajuste con 2e10, 1e11, 1e12, qgrad_ext y prediccion para qgrad (300 K)"]
    A = {}
    for er in ER:
        base = ["2e10", "1e11", "1e12", "qgrad_ext"]
        for y in ("sigma_max", "delta_sigma"):
            for xk in ("PE0", "fico_Cu0"):
                x = np.array([mean(s, er, xk) for s in base]); yy = np.array([mean(s, er, y) for s in base]); b, a = np.polyfit(x, yy, 1); res = yy - (a + b * x); se = np.sqrt((res ** 2).sum() / 2)
                pred = a + b * mean("qgrad", er, xk); A[(er, y, xk)] = (pred, mean("qgrad", er, y), se)
                L.append(f"  {er:5s} {y:11s} vs {xk:8s}: qgrad {xk}={mean('qgrad', er, xk):.5g} -> predicho {pred:.3f}, observado {mean('qgrad', er, y):.3f} (error tipico del ajuste ~{se:.3f}; diferencia {mean('qgrad', er, y)-pred:+.3f})")
    L.append("== (B) pares emparejados (Cu, CN = 12, mismo volumen propio y de vecinos): P(activo antes del pico) ico vs no-ico")
    for s in SAMPLES:
        for er in ER:
            if (s, er) in g: L.append(f"  {s:9s} {er:5s}: P(ico)={mean(s,er,'P_ico')*100:.2f}%  P(no-ico,CN12)={mean(s,er,'P_noico_CN12')*100:.2f}%  diferencia={mean(s,er,'dif')*100:+.2f} ± {sd(s,er,'dif')*100:.2f} pts  razon={mean(s,er,'razon'):.2f}  (n pares ~{mean(s,er,'n_ico_pareados'):.0f})")
    L.append("== (C) inicio de la actividad: eps donde la fraccion activa por ventana cruza 0.5 % / 1 %, y eps(sigma_max)")
    for er in ER:
        for s in SAMPLES:
            if (s, er) in g:
                o5 = [r["onset_0.5pct"] for r in g[(s, er)]]; o1 = [r["onset_1pct"] for r in g[(s, er)]]
                L.append(f"  {er:5s} {s:9s}: inicio0.5%={np.nanmean(o5):.3f} ± {np.nanstd(o5, ddof=1):.3f}  inicio1%={np.nanmean(o1):.3f} ± {np.nanstd(o1, ddof=1):.3f}  eps_pico={mean(s,er,'eps_pico'):.3f} ± {sd(s,er,'eps_pico'):.3f}")
    open(f"{M}/evidencia_ico_resumen.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
    fg, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    for s in SAMPLES:
        for er, mk in zip(ER, "s^o"):
            if (s, er) in g:
                ax[0].errorbar(mean(s, er, "fico_Cu0"), mean(s, er, "sigma_max"), yerr=sd(s, er, "sigma_max"), color=COL[s], marker=mk, mfc=COL[s] if er == "1e9" else "none", capsize=2, lw=0, elinewidth=1)
                ax[1].errorbar(mean(s, er, "PE0"), mean(s, er, "sigma_max"), yerr=sd(s, er, "sigma_max"), color=COL[s], marker=mk, mfc=COL[s] if er == "1e9" else "none", capsize=2, lw=0, elinewidth=1)
    for er in ER:
        base = ["2e10", "1e11", "1e12", "qgrad_ext"]
        for a_, xk in ((ax[0], "fico_Cu0"), (ax[1], "PE0")):
            x = np.array([mean(s_, er, xk) for s_ in base]); y = np.array([mean(s_, er, "sigma_max") for s_ in base]); b, a0 = np.polyfit(x, y, 1); xx = np.linspace(min(x.min(), mean("qgrad", er, xk)), max(x.max(), mean("qgrad", er, xk)), 10); a_.plot(xx, a0 + b * xx, color="#8a8985", lw=.8, ls="--", zorder=0)
    ax[0].set(xlabel="fracción de Cu en icosaedro (ε = 0)", ylabel="σmax (GPa)", title="σmax vs icosaedros (recta: ajuste con los 4 vidrios principales)"); ax[1].set(xlabel="Ep inicial (eV/átomo)", ylabel="σmax (GPa)", title="σmax vs energía potencial (misma convención)")
    for s in SAMPLES:
        pts = [(mean(s, er, "dif") * 100, sd(s, er, "dif") * 100) for er in ER if (s, er) in g]
        ax[2].errorbar([SAMPLES.index(s) + (i - 1) * .15 for i in range(len(pts))], [p[0] for p in pts], yerr=[p[1] for p in pts], color=COL[s], marker="o", lw=0, capsize=2)
    ax[2].axhline(0, color="#8a8985", lw=.8); ax[2].set_xticks(range(5)); ax[2].set_xticklabels(["2e10", "1e11", "1e12", "qgrad_ext", "qgrad\n(control)"]); ax[2].axvline(3.5, color="#8a8985", lw=.8, ls=":"); ax[2].set(ylabel="P(activo | ico) − P(activo | no-ico), puntos porcentuales", title="Pares emparejados (Cu, CN = 12): 1e8 / 1e9 / 1e10")
    for a in ax: a.grid(color="#e6e5e1", lw=.7); a.set_axisbelow(True); [a.spines[x].set_visible(False) for x in ("top", "right")]
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], color=COL[s], marker="o", ls="", label=(s if s != "qgrad" else "qgrad (control exploratorio, gris)")) for s in SAMPLES] + [Line2D([], [], color="#8a8985", ls="--", label="ajuste lineal con los 4 vidrios principales")] + [Line2D([], [], color="k", marker=m, ls="", mfc="none", label=f"tasa de deformación {e} s⁻¹") for m, e in (("s", "1e8"), ("^", "1e9 (relleno)"), ("o", "1e10"))]
    ax[0].legend(handles=hs, frameon=False, fontsize=8.5, loc="upper left"); ax[1].text(.02, .97, "mismos colores, marcadores y recta que a la izquierda", transform=ax[1].transAxes, fontsize=8.5, va="top")
    ax[2].text(.02, .04, "para cada vidrio, de izquierda a derecha: 1e8, 1e9, 1e10 s⁻¹\nbarras = desvío entre x, y, z del mismo vidrio", transform=ax[2].transAxes, fontsize=8.5)
    for a in ax: a.title.set_fontsize(10)
    fg.savefig(f"{M}/figs/evidencia_ico.png", dpi=135)
