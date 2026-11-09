#!/usr/bin/env python3
"""Figuras y resumen de la busqueda de STZ (lee stz_windows.csv, <caso>/analysis/stz/stz_d2min.npz, voro_frames.npz, curvas y dumps; escribe figs/stz_*.png y stz_resumen.txt).
  stz_mapas_D2min.png      mapas de D2min (rebanada de 12 A, plano x-y; traccion en x) para 1e12, 2e10 y qgrad_ext, 1e9, en 4 ventanas de d_eps = 0.01
  stz_evolucion.png        sigma(eps), fraccion de activos, tamano medio de cluster (>=2) y mayor cluster, real vs control aleatorio, por tasa de deformacion
  stz_tamanos.png          distribucion de tamanos de cluster de atomos activos (eps 0.12-0.20, 1e9) frente al control aleatorio
  stz_vs_estructura.png    actividad por atomo (ventanas activas hasta el pico) segun la estructura inicial (especie, CN, ico Cu)
Umbral: percentil 99.9 de D2min en la primera ventana (eps 0-0.01) de cada caso. Uso (anaconda + OVITO): python stz_figures.py"""
import os, csv, glob, warnings
warnings.filterwarnings("ignore")
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); F = f"{M}/figs"
SAMPLES = ["2e10", "1e11", "1e12", "qgrad_ext"]; ER = ["1e8", "1e9", "1e10"]; DIRS = "xyz"
COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad_ext": "#b5338a"}; LS = {"qgrad_ext": "--"}; GRID = "#e6e5e1"
def clean(ax):
    ax.grid(color=GRID, lw=.7); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
rel = [r for r in csv.DictReader(open(f"{M}/stz_windows.csv")) if r["umbral"].startswith("p")]
res = {r["case_id"]: r for r in csv.DictReader(open(f"{M}/resumen.csv"))}
W = defaultdict(list)
for r in rel: W[r["case_id"]].append(r)
def series(s, er, key, null=False):
    a = np.array([[float(r[("null_" + key) if null else key]) for r in sorted(W[f"{s}_300K_{er}_{d}"], key=lambda r: float(r["eps_fin"]))] for d in DIRS])
    return a.mean(0), a.std(0, ddof=1)
EPS = 0.01 * np.arange(1, 21)
def peak_idx(s, er): return int(np.clip(round(np.mean([float(res[f"{s}_300K_{er}_{d}"]["eps_sigma_max"]) for d in DIRS]) / 0.01) - 1, 0, 19))

# ---- evolucion ----
fg, axs = plt.subplots(4, 3, figsize=(16, 13), constrained_layout=True, sharex="col")
for j, er in enumerate(ER):
    for s in SAMPLES:
        cs = [f"{s}_300K_{er}_{d}" for d in DIRS]; sig = [np.loadtxt(glob.glob(f"{M}/{c}/*.curve.dat")[0], comments="#")[:199, {"x": 2, "y": 3, "z": 4}[c[-1]]] for c in cs]
        e = np.loadtxt(glob.glob(f"{M}/{cs[0]}/*.curve.dat")[0], comments="#")[:199, 1]; axs[0, j].plot(e, np.mean(sig, 0), color=COL[s], ls=LS.get(s, "-"), lw=1.5, label=s)
        for i, key in ((1, "frac_activos"), (2, "tam_medio_ge2"), (3, "mayor")):
            m, sd = series(s, er, key); axs[i, j].errorbar(EPS, m, yerr=sd, color=COL[s], ls=LS.get(s, "-"), lw=1.5, marker="o", ms=3, capsize=1.5)
            if i >= 2: axs[i, j].plot(EPS, series(s, er, key, null=True)[0], color=COL[s], ls=":", lw=1.0, alpha=.8)
        for i in range(4): axs[i, j].axvline(np.mean([float(res[c]["eps_sigma_max"]) for c in cs]), color=COL[s], lw=.8, alpha=.5)
    axs[0, j].set_title(f"tasa de deformación {er} s⁻¹", fontsize=11); axs[3, j].set_xlabel("ε (fin de la ventana de Δε = 0.01)"); axs[2, j].set_yscale("log"); axs[3, j].set_yscale("log")
for i, yl in enumerate(["σ del eje de tracción (GPa)", "fracción de átomos activos (D²min > umbral)", "tamaño medio de cluster (≥ 2 átomos)", "mayor cluster (átomos)"]): axs[i, 0].set_ylabel(yl)
for a in axs.flat: clean(a)
from matplotlib.lines import Line2D
axs[0, 0].legend(frameon=False, fontsize=9, title="vidrio")
axs[2, 0].legend(handles=[Line2D([], [], color="#444", ls="-", label="datos"), Line2D([], [], color="#444", ls=":", label="control: mismos átomos activos elegidos al azar")], frameon=False, fontsize=8.5, loc="lower right")
axs[1, 0].legend(handles=[Line2D([], [], color="#444", lw=.8, alpha=.6, label="línea vertical = ε(σmax) del vidrio (mismo color)")], frameon=False, fontsize=8.5, loc="upper left")
fg.suptitle("Actividad no afín (D²min, Δε = 0.01) y agrupamiento de átomos activos (300 K). Punteado: control con la misma cantidad de átomos elegidos al azar.\nUmbral = percentil 99.9 de la 1ª ventana de cada caso. Barras: desvío entre x, y, z. Líneas verticales: ε(σmax)", fontsize=10.5)
fg.savefig(f"{F}/stz_evolucion.png", dpi=140); plt.close(fg)

# ---- resumen numerico ----
L = []
for er in ER:
    for s in SAMPLES:
        fa = series(s, er, "frac_activos")[0]; nm = series(s, er, "tam_medio_ge2")[0]; nn = series(s, er, "tam_medio_ge2", null=True)[0]
        ma = series(s, er, "mayor")[0]; mn = series(s, er, "mayor", null=True)[0]; fc = series(s, er, "frac_en_cl_ge5")[0]; fcn = series(s, er, "frac_en_cl_ge5", null=True)[0]
        on = next((EPS[i] for i in range(1, 20) if fa[i] > 3 * max(fa[0], 1e-4)), np.nan); pk = peak_idx(s, er); st = slice(11, 20)
        L.append(f"{s:9s} {er:5s} eps_pico={0.01*(pk+1):.2f} | activos 1a ventana={fa[0]*100:.2f}% pico={fa[pk]*100:.2f}% estacionario(0.12-0.20)={fa[st].mean()*100:.2f}% | inicio(3x la 1a ventana)={on:.2f} | tam.medio>=2: pico={nm[pk]:.1f} est={nm[st].mean():.1f} (control {nn[st].mean():.1f}) | mayor est={ma[st].mean():.0f} (control {mn[st].mean():.0f}) | frac en clusters>=5: est={fc[st].mean():.2f} (control {fcn[st].mean():.2f})")
open(f"{M}/stz_resumen.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L), flush=True)

# ---- mapas y tamanos (OVITO) ----
from ovito.io import import_file
def frame_pos(c, k):
    d = import_file(glob.glob(f"{M}/{c}/*.dump")[0]).compute(k); ids = d.particles["Particle Identifier"][...]; o = np.argsort(ids); cell = np.array(d.cell[...])
    return d.particles["Position"][...][o].astype(np.float64), np.array([cell[0, 0], cell[1, 1], cell[2, 2]]), cell[:, 3]
WIN = [(0.03, "ε 0.02–0.03"), (None, "pico"), (0.10, "ε 0.09–0.10"), (0.20, "ε 0.19–0.20")]
fg, axs = plt.subplots(3, 4, figsize=(18, 12.5), constrained_layout=True)
for i, s in enumerate(["1e12", "2e10", "qgrad_ext"]):
    c = f"{s}_300K_1e9_x"; z = np.load(f"{M}/{c}/analysis/stz/stz_d2min.npz"); pk = peak_idx(s, "1e9")
    for j, (ef, lab) in enumerate(WIN):
        w = pk if ef is None else int(round(ef / 0.01)) - 1; pos, L_, org = frame_pos(c, 2 * (w + 1)); pos = (pos - org) % L_
        slab = pos[:, 2] < 12.0; a = z["d2min"][w][slab]; ax = axs[i, j]; o = np.argsort(a)
        sc = ax.scatter(pos[slab, 0][o], pos[slab, 1][o], c=np.log10(np.maximum(a[o], 0.3)), s=3.5, cmap="magma_r", vmin=0.0, vmax=np.log10(60.0), linewidths=0, rasterized=True)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f"{s} — " + (f"pico (ε 0.0{w}–{0.01*(w+1):.2f})" if ef is None else lab), fontsize=10.5, color=COL[s])
cb = fg.colorbar(sc, ax=axs, shrink=0.6, pad=0.01); cb.set_label("D²min (Å²), rebanada z < 12 Å"); cb.set_ticks(np.log10([1, 3, 10, 30, 60])); cb.set_ticklabels(["1", "3", "10", "30", "60"])
fg.suptitle("Mapas de desplazamiento no afín D²min (Δε = 0.01, 300 K, 1e9 s⁻¹; tracción en x = horizontal). Cada punto es un átomo de una rebanada de 12 Å", fontsize=11)
fg.savefig(f"{F}/stz_mapas_D2min.png", dpi=130); plt.close(fg)

from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
def csizes(pos, box, sel):
    ii = np.flatnonzero(sel); pr = cKDTree(pos[ii], boxsize=box).query_pairs(3.75, output_type="ndarray")
    return np.bincount(connected_components(coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(ii), len(ii))), directed=False)[1])
fg, axs = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True); rng = np.random.default_rng(2); extra = []
for i, s in enumerate(SAMPLES):
    c = f"{s}_300K_1e9_x"; z = np.load(f"{M}/{c}/analysis/stz/stz_d2min.npz"); allz, alln = [], []
    for w in range(11, 20):
        pos, L_, org = frame_pos(c, 2 * (w + 1)); pos = (pos - org) % L_; sel = z["d2min"][w] > float(z["thr_rel"])
        allz += list(csizes(pos, L_, sel)); alln += list(csizes(pos, L_, rng.permutation(sel)))
    for arr, ls, lab in ((allz, "-", s), (alln, ":", None)):
        v = np.sort(np.array(arr)); axs[0].loglog(v, 1 - np.arange(len(v)) / len(v), color=COL[s], ls=ls, lw=1.6, label=lab)
    v = np.array(allz); axs[1].bar(i, v[v >= 20].sum() / v.sum(), color=COL[s])
    extra.append(f"{s}: clusters (ε 0.12–0.20, x): n={len(v)}, ≥5: {(v>=5).sum()}, ≥20: {(v>=20).sum()}, ≥100: {(v>=100).sum()}, mayor {v.max()}; fracción de activos en clusters ≥20: {v[v>=20].sum()/v.sum():.2f}")
axs[0].axvspan(20, 100, color="#8a8985", alpha=.15); axs[0].text(21, 3e-4, "STZ típica\n20–100 átomos", fontsize=9)
axs[0].set(xlabel="tamaño de cluster s (átomos)", ylabel="P(tamaño ≥ s)", title="Clusters de átomos activos (ε 0.12–0.20, 1e9, x); punteado: control aleatorio"); axs[0].plot([], [], ":", color="#8a8985", label="control aleatorio (punteado)"); axs[0].legend(frameon=False)
axs[1].set(xticks=range(4), xticklabels=SAMPLES, ylabel="fracción de átomos activos en clusters ≥ 20 átomos", title="(ε 0.12–0.20)")
for a in axs: clean(a)
fg.savefig(f"{F}/stz_tamanos.png", dpi=140); plt.close(fg); print("\n".join(extra), flush=True)

# ---- actividad vs estructura inicial ----
CL = [("Cu ico", lambda t, cn, cl: (t == 2) & (cl == 1)), ("Cu no-ico", lambda t, cn, cl: (t == 2) & (cl != 1)), ("Zr CN=15", lambda t, cn, cl: (t == 1) & (cn == 15)), ("Zr otros", lambda t, cn, cl: (t == 1) & (cn != 15)), ("Ag", lambda t, cn, cl: t == 3)]
tab = {}
for s in SAMPLES:
    vals = defaultdict(list)
    for d in DIRS:
        c = f"{s}_300K_1e9_{d}"; z = np.load(f"{M}/{c}/analysis/stz/stz_d2min.npz"); v = np.load(f"{M}/{c}/analysis/structure/voro_frames.npz"); pk = peak_idx(s, "1e9")
        n_act = (z["d2min"][: pk + 2] > float(z["thr_rel"])).sum(0); t, cn, cl = v["type"], v["cn"][0], v["cls"][0]
        for nm, f in CL: vals[nm].append(n_act[f(t, cn, cl)].mean() / (pk + 2))
        for k in (11, 12, 13, 14): vals[f"cn{k}"].append(n_act[cn == k].mean() / (pk + 2))
        vals["all"].append(n_act.mean() / (pk + 2))
    tab[s] = {k: (np.mean(v), np.std(v, ddof=1)) for k, v in vals.items()}
fg, axs = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True); x = np.arange(len(CL)); wd = 0.2
for i, s in enumerate(SAMPLES): axs[0].bar(x + (i - 1.5) * wd, [tab[s][nm][0] / tab[s]["all"][0] for nm, _ in CL], wd, yerr=[tab[s][nm][1] / tab[s]["all"][0] for nm, _ in CL], color=COL[s], label=s, error_kw=dict(lw=.8, capsize=2))
axs[0].axhline(1, color="#8a8985", lw=.8); axs[0].set_xticks(x); axs[0].set_xticklabels([n for n, _ in CL]); axs[0].set_ylabel("actividad / actividad media de todos los átomos"); axs[0].set_title("Fracción de ventanas activas hasta el pico, según el entorno inicial", fontsize=10.5); axs[0].legend(frameon=False, fontsize=9)
for s in SAMPLES: axs[1].plot([11, 12, 13, 14], [tab[s][f"cn{k}"][0] / tab[s]["all"][0] for k in (11, 12, 13, 14)], color=COL[s], ls=LS.get(s, "-"), marker="o", label=s)
axs[1].legend(frameon=False, fontsize=9); axs[1].axhline(1, color="#8a8985", lw=.8); axs[1].set(xlabel="CN inicial (todos los átomos)", ylabel="actividad relativa", title="Actividad hasta el pico vs CN inicial"); axs[1].set_xticks([11, 12, 13, 14])
for i, s in enumerate(SAMPLES): axs[2].bar(i, tab[s]["all"][0] * 100, color=COL[s], yerr=tab[s]["all"][1] * 100, error_kw=dict(lw=.8, capsize=2))
axs[2].set(xticks=range(4), xticklabels=SAMPLES, ylabel="% de ventanas activas por átomo (hasta el pico)", title="Actividad media hasta el pico (1e9)")
for a in axs: clean(a)
fg.suptitle("¿Dónde nacen las STZ? Actividad no afín acumulada hasta el pico según la estructura inicial (media de x, y, z; barras = desvío entre direcciones)", fontsize=11)
fg.savefig(f"{F}/stz_vs_estructura.png", dpi=140); plt.close(fg)
T = "\n".join(f"{s}: " + "  ".join(f"{nm}={tab[s][nm][0]/tab[s]['all'][0]:.2f}" for nm, _ in CL) + "  | " + " ".join(f"CN{k}={tab[s][f'cn{k}'][0]/tab[s]['all'][0]:.2f}" for k in (11, 12, 13, 14)) + f" | media={tab[s]['all'][0]*100:.2f}%" for s in SAMPLES)
open(f"{M}/stz_resumen.txt", "a").write("\n== tamanos de cluster (eps 0.12-0.20, 1e9, x)\n" + "\n".join(extra) + "\n\n== actividad relativa por entorno inicial (hasta el pico, 1e9)\n" + T + "\n"); print(T)
