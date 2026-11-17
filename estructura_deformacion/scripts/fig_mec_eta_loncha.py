"""Bloque 8 — ¿la deformación atómica sigue a la estructura inicial? Mapas de η (von Mises, acumulada desde ε = 0) en la
MISMA loncha que los mapas de q6 (z < 12 Å), tracción en x a 300 K y 1e9 s^-1, y una prueba cuantitativa por celdas.
Muestras: 1e12 (unificada, cuando esté), 1e11, 2e10, qgrad_ext (sin qgrad).

Lee  : mechanical/matrix/<caso>/analysis/eta_loncha.npz (scripts/eta_loncha.py: η a ε = 0,06 y 0,10; posiciones del cuadro 0)
       qgrad/laves/features/<muestra>/17.npz (índice de Voronoi y vecinos del vidrio a 300 K; mismos ids)
Sale : estructura_deformacion/figs/mec_eta_loncha_<muestra>_<eps>.png   mapa de η, contorno negro = Cu-ico en clusters de 5 o más (como est_mapa_q6)
       estructura_deformacion/figs/mec_eta_vs_ico_<eps>.png              η media de cada celda (~10,5 Å) / η media global, contra la fracción de
                                                           Cu icosaédrico de la celda (por quintiles de esa fracción)
Chequeo: Spearman átomo a átomo η–q6 y celda a celda η–fracción ico; η de los Cu-ico en clusters de 5 o más / η de los otros Cu.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr

RAIZ = Path("/scratch/negishi/aaabt/bmg_cuzrag")
DEST = RAIZ / "estructura_deformacion" / "figs"   # análisis en curso: NO va a la bitácora (pedido 17/11/2026)
sys.path.insert(0, str(RAIZ / "estilo"))
import bmg_estilo as be  # noqa: E402
import matplotlib as mpl  # noqa: E402

M, SA = "mechanical/matrix", "structural_analysis"
MUESTRAS = {"1e12": ("1e12u_300K_1e9_x", f"{SA}/qgrad/laves/features/1e12u/17.npz"),
            "1e11": ("1e11_300K_1e9_x", f"{SA}/qgrad/laves/features/1e11/17.npz"),
            "2e10": ("2e10_300K_1e9_x", f"{SA}/qgrad/laves/features/2e10/17.npz"),
            "qgrad_ext": ("qgrad_ext_300K_1e9_x", "qgrad_ext/laves/features/qgrad_ext/17r.npz")}
EPS = {"006": (0.06, 0.25), "010": (0.10, 0.40)}      # ε y máximo de la escala de color
SLAB, NMIN, NV = 12.0, 5, 10                            # loncha, cluster mínimo, celdas por lado (10^3 celdas de ~10,5 Å)
ICO = np.array([0, 0, 12, 0])

curvas = {e: {} for e in EPS}
for k, (caso, pf) in MUESTRAS.items():
    f = RAIZ / M / caso / "analysis" / "eta_loncha.npz"
    if not f.exists():
        print(f"{k}: falta {f.relative_to(RAIZ)} (tracción todavía no corrida) -> salteado")
        continue
    E = dict(np.load(f)); F = dict(np.load(RAIZ / pf))
    t = E["type"]; assert (F["type"] == t).all(), k
    n = len(t); x = E["x0"] % E["L0"]; L = E["L0"]
    cu = t == 2; ico = cu & (F["vindex"] == ICO).all(1)
    src = np.repeat(np.arange(n), np.diff(F["neigh_ptr"])); idx = F["neigh_idx"]; m = ico[src] & ico[idx]
    _, cl = connected_components(csr_matrix((np.ones(m.sum()), (src[m], idx[m])), shape=(n, n)), directed=False)
    grande = ico & (np.bincount(cl, weights=ico.astype(float))[cl] >= NMIN)
    vox = (np.floor(x / L * NV).astype(int) % NV) @ np.array([NV * NV, NV, 1])
    cnt_cu = np.bincount(vox, weights=cu.astype(float), minlength=NV ** 3)
    fico = np.bincount(vox, weights=ico.astype(float), minlength=NV ** 3) / np.maximum(cnt_cu, 1)
    s = x[:, 2] < SLAB
    for e, (eps, vmax) in EPS.items():
        eta = E[f"eta_{e}"]
        eta_v = np.bincount(vox, weights=eta, minlength=NV ** 3) / np.bincount(vox, minlength=NV ** 3)
        rho_at = spearmanr(eta, F["q6_12"])[0]; rho_v = spearmanr(fico, eta_v)[0]
        q = np.quantile(fico, np.linspace(0, 1, 6)); b = np.clip(np.searchsorted(q, fico, side="right") - 1, 0, 4)
        curvas[e][k] = (np.array([fico[b == i].mean() for i in range(5)]), np.array([eta_v[b == i].mean() for i in range(5)]) / eta.mean(),
                        np.array([eta_v[b == i].std() / np.sqrt((b == i).sum()) for i in range(5)]) / eta.mean())
        print(f"{k:9s} ε={eps:.2f}: η media {eta.mean():.4f}; Spearman átomo η–q6 = {rho_at:+.3f}; celdas η–f_ico = {rho_v:+.3f}; "
              f"η(Cu-ico en clusters ≥{NMIN}) / η(otros Cu) = {eta[grande].mean() / eta[cu & ~grande].mean():.3f}")
        norm = mpl.colors.Normalize(0, vmax)
        o = np.argsort(eta[s])
        fig, ax = be.figura("cuadrada")
        sc = ax.scatter(x[s, 0][o], x[s, 1][o], c=eta[s][o], s=0.5 + 6 * np.clip(eta[s][o] / vmax, 0, 1) ** 2,
                        cmap=be.CMAP_SEC, norm=norm, linewidths=0)
        g = s & grande
        ax.scatter(x[g, 0], x[g, 1], s=12, facecolors="none", edgecolors="k", linewidths=0.7)   # negro: el blanco no contrastaba (17/11)
        ax.set_aspect("equal"); ax.set_xlim(0, L[0]); ax.set_ylim(0, L[1]); ax.set_xticks([]); ax.set_yticks([])
        cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.046); cb.set_label(r"$\eta$")
        be.guardar(fig, f"mec_eta_loncha_{k}_{e}", script=__file__, destino=DEST, fuentes=[f"{M}/{caso}/analysis/eta_loncha.npz", pf])

for e, (eps, _) in EPS.items():
    fig, ax = be.figura("media")
    ax.axhline(1, color="0.6", lw=0.6, zorder=0)
    for k, (fx, fy, fe) in curvas[e].items():
        ax.errorbar(100 * fx, fy, yerr=fe, ms=4, lw=1.0, capsize=0, **be.linea(k))
    ax.set_xlabel(r"Cu $\langle0,0,12,0\rangle$ en la celda (% de los Cu)")
    ax.set_ylabel(r"$\eta$ de la celda / $\eta$ media")
    be.leyenda(ax, loc="upper right", fontsize=7)
    be.guardar(fig, f"mec_eta_vs_ico_{e}", script=__file__, destino=DEST, fuentes=[f"{M}/*/analysis/eta_loncha.npz"])
