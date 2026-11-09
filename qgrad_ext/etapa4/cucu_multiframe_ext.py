#!/usr/bin/env python3
"""g_Cu-Cu(r) multi-frame a 300 K para las 4 curvas, sobre los 10 frames (1 por ~ns) de las relajaciones NPT de 10 ns (dumps
relax300K_*.dump; se descarta el frame inicial). Objetivo: ver si la profundidad del valle Cu-Cu de qgrad a 300 K (0.400, 1 snapshot)
vs 2e10 (0.353) es senal o fluctuacion. Misma grilla que LAMMPS (500 bins, 0-10 A), mismas ventanas que cucu_split_4curvas.py.
g(r) por conteo de pares periodico (cKDTree); se VALIDA contra el g_CuCu de LAMMPS (rdf_sametype_all.csv, run 0 del restart final).
Salidas: cucu_multiframe_ext.csv, ../fig_R_CuCu_multiframe_4curvas.png.  Solo lectura de los dumps.
COPIA de structural_analysis/qgrad/etapa4/rdf_sametype/cucu_multiframe.py con qgrad_ext agregado (su "1 snapshot" sale de snaps/17r.dump con el mismo conteo)."""
import csv
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent; B = Path("/scratch/negishi/aaabt/bmg_cuzrag"); REFCSV = B / "structural_analysis/qgrad/etapa4/rdf_sametype/rdf_sametype_all.csv"
DUMPS = {"2e10": B / "rate_2e10/prod_2e10/relax300K_N65536_2e10.dump", "1e11": B / "rate_1e11/prod_1e11/relax300K_N65536_1e11.dump",
         "1e12": B / "rate_1e12_reference/relax300K_N65536_1e12.dump", "qgrad": B / "structural_analysis/qgrad/relax300K/relax300K_N65536_qgrad.dump",
         "qgrad_ext": B / "qgrad_ext/relax300K_N65536_gradext.dump"}
W1, WD, W2 = (4.0, 4.9), (4.7, 5.3), (5.0, 5.7)
EDGES = np.linspace(0, 10, 501); RC = 0.5 * (EDGES[1:] + EDGES[:-1])
SHELL = 4 / 3 * np.pi * (EDGES[1:] ** 3 - EDGES[:-1] ** 3)
def frames(path):
    txt = open(path).read().split("ITEM: TIMESTEP")[1:]
    for t in txt:
        L = t.strip().split("\n"); n = int(L[2]); b = np.array([l.split()[:2] for l in L[4:7]], float)
        a = np.array([l.split()[:5] for l in L[8:8 + n]], float); a = a[np.argsort(a[:, 0])]
        yield int(L[0]), a[:, 1].astype(int), a[:, 2:5] - b[:, 0], (b[:, 1] - b[:, 0])[0]
def gcucu(t, x, L):
    cu = x[t == 2] % L; n = len(cu); tr = cKDTree(cu, boxsize=L)
    pr = tr.query_pairs(10.0, output_type="ndarray")                # pares (i<j) a < 10 A, periodico
    dv = cu[pr[:, 0]] - cu[pr[:, 1]]; d = np.linalg.norm(dv - L * np.round(dv / L), axis=1)
    h = np.histogram(d, bins=EDGES)[0] * 2.0                       # pares ordenados
    return h / (n * (n / L ** 3) * SHELL)
def depth(g):
    def ext(w, f): m = (RC > w[0]) & (RC < w[1]); i = f(g[m]); return RC[m][i], g[m][i]
    r1, g1 = ext(W1, np.argmax); rd, gd = ext(WD, np.argmin); r2, g2 = ext(W2, np.argmax)
    return (g1 + g2) / 2 - gd, r2
ref = {}
for r in csv.DictReader(open(REFCSV)):
    if int(float(r["T"])) == 300: ref.setdefault(r["rate"], []).append(float(r["g_CuCu"]))
def snap_g(label):
    for step, t, x, L in frames(str(HERE.parent / "snaps" / f"{label}.dump")): return gcucu(t, x, L)
ref["qgrad_ext"] = list(snap_g("17r"))
out, G, GF = [], {}, {}
for k, p in DUMPS.items():
    fr = list(frames(p)); assert len(fr) == 11, (k, len(fr))
    gs = []
    for i, (step, t, x, L) in enumerate(fr):
        g = gcucu(t, x, L); gs.append(g)
        if i >= 1: d, r2 = depth(g); out.append([k, i, step, d, r2])
    # validacion: ultimo frame == restart final -> g de LAMMPS
    if k in ref and len(ref[k]) == 500:
        print(f"validacion {k:6s}: max|g_mio(ultimo frame) - g_LAMMPS(restart final)| = {np.max(np.abs(gs[-1] - np.array(ref[k]))):.4f}   (depth mio {depth(gs[-1])[0]:.4f} vs LAMMPS {depth(np.array(ref[k]))[0]:.4f})")
    G[k] = np.mean(gs[1:], axis=0); GF[k] = np.array(gs[1:])
with open(HERE / "cucu_multiframe_ext.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "frame", "step", "depth", "r_peak2"]); w.writerows([[a, b, c, f"{d:.4f}", f"{e:.2f}"] for a, b, c, d, e in out])
print(f"\n{'curva':6s} {'depth: media de 10 frames':>26s} {'std':>7s} {'error std (std/sqrt10)':>24s} {'depth del g promediado':>24s} {'1 snapshot (restart)':>21s}")
S = {}
for k in DUMPS:
    d = np.array([o[3] for o in out if o[0] == k]); dm = depth(G[k])[0]
    one = depth(np.array(ref[k]))[0] if k in ref else float("nan"); S[k] = (d.mean(), d.std(ddof=1), dm)
    print(f"{k:6s} {d.mean():26.4f} {d.std(ddof=1):7.4f} {d.std(ddof=1)/np.sqrt(len(d)):24.4f} {dm:24.4f} {one:21.4f}")
print("\nDiferencia de profundidad respecto de 2e10 (media de 10 frames; error = sqrt(sem_a^2+sem_b^2), SOLO ruido termico entre frames):")
for k in ("1e11", "1e12", "qgrad", "qgrad_ext"):
    dd = S[k][0] - S["2e10"][0]; e = np.hypot(S[k][1], S["2e10"][1]) / np.sqrt(10)
    print(f"  {k:6s} - 2e10 = {dd:+.4f} +- {e:.4f}  ({dd/e:+.1f} sigma)")
# depth del g PROMEDIADO (menos sesgado que el de un frame: el maximo de picos ruidosos sesga hacia arriba) con bootstrap sobre los 10 frames
rng = np.random.default_rng(0); BS = {}
for k in DUMPS:
    BS[k] = np.array([depth(GF[k][rng.integers(0, 10, 10)].mean(0))[0] for _ in range(2000)])
print("\nDepth del g promediado (10 frames), media y error bootstrap (SOLO ruido termico entre frames):")
for k in DUMPS: print(f"  {k:6s} {depth(G[k])[0]:.4f} +- {BS[k].std():.4f}")
print("Diferencias respecto de 2e10 (bootstrap pareado por indice de remuestreo no aplica: curvas independientes):")
for k in ("1e11", "1e12", "qgrad", "qgrad_ext"):
    dd = depth(G[k])[0] - depth(G["2e10"])[0]; e = np.hypot(BS[k].std(), BS["2e10"].std())
    print(f"  {k:6s} - 2e10 = {dd:+.4f} +- {e:.4f}  ({dd/e:+.1f} sigma)")
with open(HERE / "cucu_multiframe_resumen_ext.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "depth_media_frames", "std_frames", "depth_g_promediado", "err_bootstrap", "depth_1_snapshot_restart"])
    for k in DUMPS: w.writerow([k, f"{S[k][0]:.4f}", f"{S[k][1]:.4f}", f"{depth(G[k])[0]:.4f}", f"{BS[k].std():.4f}", f"{depth(np.array(ref[k]))[0]:.4f}"])
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2)); COL = {"1e12": "#d62728", "1e11": "#e6b800", "2e10": "#1f3fd6", "qgrad": "#7b3fbf", "qgrad_ext": "#0b0b0b"}
for j, k in enumerate(DUMPS):
    d = [o[3] for o in out if o[0] == k]; ax[0].plot(np.full(len(d), j) + np.linspace(-.15, .15, len(d)), d, "o", ms=4, color=COL[k], alpha=.6)
    ax[0].errorbar(j + .3, S[k][0], yerr=S[k][1], fmt="s", color=COL[k], capsize=4, ms=7)
ax[0].set_xticks(range(5)); ax[0].set_xticklabels(list(DUMPS)); ax[0].set(ylabel="profundidad del valle Cu-Cu", title="10 frames por curva (puntos), media +- std (cuadrado)")
for k in DUMPS: ax[1].plot(RC, G[k], color=COL[k], lw=1.3, label=k)
ax[1].set(xlim=(3.8, 6.4), xlabel="r (A)", ylabel="g_Cu-Cu(r) promedio de 10 frames", title="g_Cu-Cu 300 K, promedio multi-frame"); ax[1].legend(fontsize=8)
for a in ax: a.grid(alpha=.25)
fig.tight_layout(); fig.savefig(HERE / "fig_R_CuCu_multiframe_5curvas.png", dpi=160)
