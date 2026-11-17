#!/usr/bin/env python3
"""Evolución de las STZ y de la relación deformación–estructura con la deformación (17/11/2026). Tracción en x, 300 K, 1e9 s^-1.
Muestras: 1e12 (unificada; se agrega sola cuando exista su eta_loncha.npz), 1e11, 2e10, qgrad_ext (sin qgrad).

Lee (solo lectura):
  mechanical/matrix/<caso>/analysis/eta_loncha.npz   (eta_loncha.py: η acumulada a ε = 0,02 ... 0,20; posiciones de ε = 0)
  mechanical/matrix/<caso>/analysis/stz/stz_d2min.npz (stz_analysis.py: D²min por ventana de Δε = 0,01, referencia = inicio de
                                                       la ventana; umbral thr_rel = p99,9 de la 1.ª ventana del caso)
  structural_analysis/qgrad/laves/features/<muestra>/17.npz (índice y vecinos de Voronoi del vidrio a 300 K; mismos ids)
Salidas (estructura_deformacion/):
  figs/stz_loncha_<muestra>_eps<ε>.png  activos de la ventana que termina en ε (rojo) sobre la loncha z < 12 Å (gris), contorno negro =
                                        Cu-ico en clusters de 5 o más; ε = 0,03 0,05 0,07 0,10 0,15 0,20
  figs/eta_loncha_<muestra>_eps<ε>.png  η acumulada hasta ε (ε = 0,02 ... 0,20), mismo contorno
  figs/evol_*.png                       curvas contra ε (ver abajo); figs/hoja_*.png: hojas de contacto solo para mirar
  datos/evolucion_deformacion.csv       una fila por muestra y ε
Uso: PYTHONDONTWRITEBYTECODE=1 /scratch/negishi/aaabt/venv_structural/bin/python evolucion_deformacion.py
"""
import csv, sys
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr

R = Path("/scratch/negishi/aaabt/bmg_cuzrag"); M = R / "mechanical/matrix"; ED = R / "estructura_deformacion"
sys.path.insert(0, str(R / "estilo"))
import bmg_estilo as be  # noqa: E402
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

SA = R / "structural_analysis"
MUESTRAS = {"1e12": ("1e12u_300K_1e9_x", SA / "qgrad/laves/features/1e12u/17.npz"),
            "1e11": ("1e11_300K_1e9_x", SA / "qgrad/laves/features/1e11/17.npz"),
            "2e10": ("2e10_300K_1e9_x", SA / "qgrad/laves/features/2e10/17.npz"),
            "qgrad_ext": ("qgrad_ext_300K_1e9_x", R / "qgrad_ext/laves/features/qgrad_ext/17r.npz")}
EPS_ETA = [0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20]
EPS_STZ_MAPA = [0.03, 0.05, 0.07, 0.10, 0.15, 0.20]
SLAB, NMIN, NV, ICO = 12.0, 5, 10, np.array([0, 0, 12, 0])
FIGS = ED / "figs"


def clusters(mask, src, idx, n):
    m = mask[src] & mask[idx]
    _, cl = connected_components(csr_matrix((np.ones(m.sum()), (src[m], idx[m])), shape=(n, n)), directed=False)
    return cl, np.bincount(cl, weights=mask.astype(float))


def mapa(nombre, x, s, color, tam, cmap, norm, contorno, etq, fuentes, gris=None):
    fig, ax = be.figura("cuadrada")
    if gris is not None:
        ax.scatter(x[s & ~gris, 0], x[s & ~gris, 1], s=0.4, color="0.8", linewidths=0)
        sel = s & gris
    else:
        sel = s
    o = np.argsort(color[sel])
    sc = ax.scatter(x[sel, 0][o], x[sel, 1][o], c=color[sel][o], s=tam[sel][o], cmap=cmap, norm=norm, linewidths=0)
    g = s & contorno
    ax.scatter(x[g, 0], x[g, 1], s=12, facecolors="none", edgecolors="k", linewidths=0.7)
    ax.set_aspect("equal"); ax.set_xlim(0, L[0]); ax.set_ylim(0, L[1]); ax.set_xticks([]); ax.set_yticks([])
    cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.046); cb.set_label(etq)
    be.guardar(fig, nombre, script=__file__, destino=FIGS, fuentes=fuentes)
    plt.close(fig)


filas, disponibles = [], []
for k, (caso, pf) in MUESTRAS.items():
    fe, fd = M / caso / "analysis/eta_loncha.npz", M / caso / "analysis/stz/stz_d2min.npz"
    if not (fe.exists() and fd.exists()):
        print(f"{k}: faltan datos ({caso}) -> salteado"); continue
    disponibles.append(k)
    E = dict(np.load(fe)); Z = dict(np.load(fd)); F = dict(np.load(pf))
    t = E["type"]; n = len(t); assert (F["type"] == t).all() and (np.sort(Z["id"]) == E["id"]).all()
    d2 = Z["d2min"][:, np.argsort(Z["id"])]; thr = float(Z["thr_rel"]); efin = Z["eps_fin"]
    L = E["L0"]; x = E["x0"] % L; s = x[:, 2] < SLAB
    cu = t == 2; ico = cu & (F["vindex"] == ICO).all(1)
    src = np.repeat(np.arange(n), np.diff(F["neigh_ptr"])); idx = F["neigh_idx"]
    cl, tam_cl = clusters(ico, src, idx, n); grande = ico & (tam_cl[cl] >= NMIN)
    vox = (np.floor(x / L * NV).astype(int) % NV) @ np.array([NV * NV, NV, 1])
    cnt = np.bincount(vox, minlength=NV ** 3); fico = np.bincount(vox, weights=ico.astype(float), minlength=NV ** 3) / np.maximum(
        np.bincount(vox, weights=cu.astype(float), minlength=NV ** 3), 1)
    fuentes = [str(fe.relative_to(R)), str(fd.relative_to(R)), str(Path(pf).relative_to(R))]
    # ---- η acumulada por nivel
    for e in EPS_ETA:
        eta = E[f"eta_{round(1000 * e):03d}"]
        eta_v = np.bincount(vox, weights=eta, minlength=NV ** 3) / cnt
        r = dict(muestra=k, eps=e, tipo="eta", valor_medio=eta.mean(),
                 eta_ico_cl_sobre_cu=eta[grande].mean() / eta[cu & ~grande].mean(),
                 eta_ico_sobre_cu=eta[ico].mean() / eta[cu & ~ico].mean(),
                 spearman_celdas=spearmanr(fico, eta_v)[0])
        filas.append(r)
        vmax = 0.6 if e >= 0.15 else 0.4 if e >= 0.08 else 0.2
        mapa(f"eta_loncha_{k}_eps{round(1000 * e):03d}", x, s, eta, 0.5 + 6 * np.clip(eta / vmax, 0, 1) ** 2, be.CMAP_SEC,
             mpl.colors.Normalize(0, vmax), grande, r"$\eta$", fuentes)
    # ---- STZ por ventana de Δε = 0,01
    for w in range(len(efin)):
        act = d2[w] > thr
        cla, tam_a = clusters(act, src, idx, n); tams = tam_a[np.unique(cla[act])] if act.any() else np.array([0])
        r = dict(muestra=k, eps=round(float(efin[w]), 3), tipo="stz", valor_medio=act.mean(),
                 p_act_ico_sobre_cu=act[ico].mean() / max(act[cu & ~ico].mean(), 1e-9),
                 p_act_ico_cl_sobre_cu=act[grande].mean() / max(act[cu & ~grande].mean(), 1e-9),
                 stz_tam_medio_ge2=tams[tams >= 2].mean() if (tams >= 2).any() else 1.0, stz_mayor=int(tams.max()),
                 spearman_celdas=spearmanr(fico, np.bincount(vox, weights=act.astype(float), minlength=NV ** 3) / cnt)[0])
        filas.append(r)
        if round(float(efin[w]), 2) in EPS_STZ_MAPA:
            col = np.log10(np.clip(d2[w], thr, 10 * thr) / thr)
            mapa(f"stz_loncha_{k}_eps{round(1000 * efin[w]):03d}", x, s, col, np.full(n, 3.0), mpl.colormaps["Reds"].resampled(256),
                 mpl.colors.Normalize(-0.3, 1), grande, r"$\log_{10}(D^2_{\min}/D^2_{\mathrm{umbral}})$", fuentes, gris=act)
    print(f"{k}: listo")

cols = ["muestra", "eps", "tipo", "valor_medio", "eta_ico_cl_sobre_cu", "eta_ico_sobre_cu", "p_act_ico_sobre_cu",
        "p_act_ico_cl_sobre_cu", "stz_tam_medio_ge2", "stz_mayor", "spearman_celdas"]
with open(ED / "datos/evolucion_deformacion.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in filas:
        w.writerow({c: (f"{r[c]:.4g}" if isinstance(r.get(c), float) else r.get(c, "")) for c in cols})

# ---- curvas contra ε
CURVAS = [("evol_stz_activos", "stz", "valor_medio", 100, "átomos activos en la ventana (%)", None),
          ("evol_stz_tamano", "stz", "stz_tam_medio_ge2", 1, "tamaño medio de cluster STZ (átomos)", None),
          ("evol_stz_ico", "stz", "p_act_ico_sobre_cu", 1, "P(activo | Cu-ico) / P(activo | otro Cu)", 1.0),
          ("evol_eta_ico", "eta", "eta_ico_cl_sobre_cu", 1, r"$\eta$(Cu-ico en clusters) / $\eta$(otros Cu)", 1.0),
          ("evol_eta_celdas", "eta", "spearman_celdas", 1, r"Spearman celdas: $\eta$ contra fracción ico", 0.0)]
for nombre, tipo, col, esc, ylab, ref in CURVAS:
    fig, ax = be.figura("media")
    if ref is not None:
        ax.axhline(ref, color="0.6", lw=0.6, zorder=0)
    ax.axvspan(0.055, 0.068, color="0.9", lw=0, zorder=0)
    for k in disponibles:
        rr = [r for r in filas if r["muestra"] == k and r["tipo"] == tipo]
        ax.plot([r["eps"] for r in rr], [esc * r[col] for r in rr], ms=3, lw=1.0, **be.linea(k))
    ax.set_xlabel(r"$\varepsilon$"); ax.set_ylabel(ylab); ax.set_xlim(0, 0.205)
    be.leyenda(ax, loc="best", fontsize=7)
    be.guardar(fig, nombre, script=__file__, destino=FIGS, fuentes=["estructura_deformacion/datos/evolucion_deformacion.csv"])
    plt.close(fig)

# ---- hojas de contacto (solo para mirar; no son figuras para el informe)
for pref, eps in (("stz_loncha", EPS_STZ_MAPA), ("eta_loncha", EPS_ETA)):
    fig, axs = plt.subplots(len(disponibles), len(eps), figsize=(2.2 * len(eps), 2.0 * len(disponibles)), squeeze=False)
    for i, k in enumerate(disponibles):
        for j, e in enumerate(eps):
            a = axs[i, j]; a.imshow(plt.imread(FIGS / f"{pref}_{k}_eps{round(1000 * e):03d}.png")); a.axis("off")
            if i == 0: a.set_title(f"ε = {e:.2f}".replace(".", ","), fontsize=9)
        axs[i, 0].text(-0.05, 0.5, k, transform=axs[i, 0].transAxes, rotation=90, va="center", ha="right", fontsize=9)
    fig.tight_layout(); fig.savefig(FIGS / f"hoja_{pref}.png", dpi=150); plt.close(fig)
print("hojas de contacto listas")
