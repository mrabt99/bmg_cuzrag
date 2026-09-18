#!/usr/bin/env python3
"""Tarea F, punto 4 -- mapa (q4, q6) por nnn, equivalente a la Fig. 8 del
paper, ANTES de proponer ninguna regla de corte para Z14 / Z15.

La pregunta que tiene que contestar el mapa es una sola: cerca del punto
ideal de Z14 (nnn=14) / Z15 (nnn=15), hay una cola o un agrupamiento
separable del grueso del vidrio, o el ideal cae adentro de la nube?

Diagnosticos por nnn:
  * densidad 2D (log) en (q4, q6), un panel por especie central;
  * contornos de las 3 especies superpuestos;
  * perfil radial de cuentas alrededor del ideal: si hubiera un
    agrupamiento propio, el perfil tendria un bump a radio chico; si el
    ideal es solo el borde/interior de la nube, decae monotono;
  * marginal de q4 dentro de una ventana de q6 alrededor del ideal
    (misma construccion que se uso para Z16), para ver si un umbral 1D
    tipo el de icosaedrico tendria de donde agarrarse.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter

from boo_types import (SPECIES, SPECIES_COLOR, IDEAL, load_boo_dump, features)

DUMP = "../boo/2e10/17.boo.dump"
CASES = [(14, "Z14"), (15, "Z15")]
Q6_HALFWIDTH = 0.01          # misma ventana que usa el criterio de Z16

data = load_boo_dump(DUMP)
types = data[:, 1].astype(int)

def ramp(color):
    return LinearSegmentedColormap.from_list("r", ["#ffffff", color])

for nnn, tname in CASES:
    ref = IDEAL[tname]
    q4, q6, _w4, w6 = features(data, nnn)

    print(f"\n{'='*66}")
    print(f"=== nnn={nnn}, ideal {tname}: q4={ref['q4']}, q6={ref['q6']}, "
          f"What6={ref['w6']} ===")
    print(f"{'especie':8s} {'q4 p5':>8s} {'q4 p50':>8s} {'q4 p95':>8s} "
          f"{'q6 p5':>8s} {'q6 p50':>8s} {'q6 p95':>8s}")
    for t, name in SPECIES.items():
        m = types == t
        print(f"{name:8s} " + " ".join(f"{v:8.4f}" for v in (
            np.percentile(q4[m], 5), np.median(q4[m]), np.percentile(q4[m], 95),
            np.percentile(q6[m], 5), np.median(q6[m]), np.percentile(q6[m], 95))))

    # donde cae el ideal dentro de la distribucion de cada especie
    print(f"\npercentil en el que cae el ideal {tname} dentro de cada especie:")
    print(f"{'especie':8s} {'pct(q4<=ideal)':>16s} {'pct(q6<=ideal)':>16s}")
    for t, name in SPECIES.items():
        m = types == t
        print(f"{name:8s} {100.0*(q4[m] <= ref['q4']).mean():16.3f} "
              f"{100.0*(q6[m] <= ref['q6']).mean():16.3f}")

    # --- perfil radial alrededor del ideal en (q4, q6) -----------------
    r = np.sqrt((q4 - ref["q4"]) ** 2 + (q6 - ref["q6"]) ** 2)
    print(f"\nperfil radial en (q4,q6) alrededor del ideal {tname} "
          f"-- densidad por unidad de area (cuentas / 2*pi*r*dr):")
    edges = np.linspace(0, 0.30, 31)
    print(f"{'r':>12s} {'Zr':>10s} {'Cu':>10s} {'Ag':>10s}")
    for lo, hi in zip(edges[:-1], edges[1:]):
        area = np.pi * (hi ** 2 - lo ** 2)
        vals = []
        for t in SPECIES:
            m = types == t
            vals.append(((r[m] >= lo) & (r[m] < hi)).sum() / area)
        if hi <= 0.16:
            print(f"{lo:5.3f}-{hi:5.3f} " + " ".join(f"{v:10.1f}" for v in vals))

    # --- marginal de q4 dentro de la ventana de q6 ---------------------
    win = np.abs(q6 - ref["q6"]) <= Q6_HALFWIDTH
    print(f"\nventana |q6 - {ref['q6']}| <= {Q6_HALFWIDTH}: {win.sum()} atomos "
          f"({100.0*win.mean():.2f} % del sistema)")
    for t, name in SPECIES.items():
        m = win & (types == t)
        print(f"    {name:2s}: {m.sum():5d} ({100.0*m.sum()/win.sum():5.2f} % "
              f"del subconjunto)")
    q4w = q4[win]
    print(f"    q4 dentro de la ventana: min={q4w.min():.4f} "
          f"p5={np.percentile(q4w,5):.4f} p50={np.median(q4w):.4f} "
          f"p95={np.percentile(q4w,95):.4f}   (ideal q4={ref['q4']})")

    # --- figura -------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.4))
    lo4, hi4 = 0.0, np.percentile(q4, 99.9)
    lo6, hi6 = np.percentile(q6, 0.1), np.percentile(q6, 99.9)

    for k, (t, name) in enumerate(SPECIES.items()):
        ax = axes[0, k]
        m = types == t
        hb = ax.hexbin(q4[m], q6[m], gridsize=85, bins="log",
                       extent=(lo4, hi4, lo6, hi6), cmap=ramp(SPECIES_COLOR[name]),
                       linewidths=0)
        ax.plot(ref["q4"], ref["q6"], "*", ms=18, color="#e34948",
                markeredgecolor="#0b0b0b", markeredgewidth=0.7, zorder=5)
        ax.set_xlim(lo4, hi4); ax.set_ylim(lo6, hi6)
        ax.set_xlabel(f"q4 (nnn={nnn})"); ax.set_ylabel(f"q6 (nnn={nnn})")
        ax.set_title(f"{name}-centrado (n={m.sum()})")
        fig.colorbar(hb, ax=ax, label="cuentas (log)")

    ax = axes[1, 0]
    for t, name in SPECIES.items():
        m = types == t
        H, xe, ye = np.histogram2d(q4[m], q6[m], bins=70,
                                   range=[[lo4, hi4], [lo6, hi6]])
        # suavizado gaussiano antes de contornear: sin esto el nivel del 2%
        # persigue ruido de binning y los contornos salen como confeti
        H = gaussian_filter(H, sigma=1.8)
        H = H / H.max()
        ax.contour(0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:]), H.T,
                   levels=[0.02, 0.15, 0.5], colors=SPECIES_COLOR[name],
                   linewidths=1.3)
        ax.plot([], [], color=SPECIES_COLOR[name], lw=1.6, label=name)
    ax.plot(ref["q4"], ref["q6"], "*", ms=18, color="#e34948",
            markeredgecolor="#0b0b0b", markeredgewidth=0.7, zorder=5,
            label=f"ideal {tname}")
    ax.set_xlim(lo4, hi4); ax.set_ylim(lo6, hi6)
    ax.set_xlabel(f"q4 (nnn={nnn})"); ax.set_ylabel(f"q6 (nnn={nnn})")
    ax.set_title("contornos por especie (2%, 15%, 50% del pico)")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.15)

    ax = axes[1, 1]
    centers = 0.5 * (edges[:-1] + edges[1:])
    for t, name in SPECIES.items():
        m = types == t
        dens = [((r[m] >= lo) & (r[m] < hi)).sum() / (np.pi * (hi ** 2 - lo ** 2))
                for lo, hi in zip(edges[:-1], edges[1:])]
        ax.plot(centers, dens, lw=1.8, color=SPECIES_COLOR[name], label=name)
    ax.set_xlabel(f"distancia radial al ideal {tname} en (q4, q6)")
    ax.set_ylabel("densidad por unidad de area")
    ax.set_title("perfil radial: un bump a r chico = agrupamiento propio")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.15)

    ax = axes[1, 2]
    bb = np.linspace(0, max(np.percentile(q4w, 99.5), ref["q4"] * 1.5), 70)
    for t, name in SPECIES.items():
        m = win & (types == t)
        if m.sum() == 0:
            continue
        ax.hist(q4[m], bins=bb, histtype="step", lw=1.8,
                color=SPECIES_COLOR[name], label=f"{name} (n={m.sum()})")
    ax.axvline(ref["q4"], color="#e34948", ls="--", lw=1.4,
               label=f"ideal q4={ref['q4']}")
    ax.set_xlabel(f"q4 (nnn={nnn})")
    ax.set_ylabel("cuentas")
    ax.set_title(f"q4 dentro de |q6-{ref['q6']}|<={Q6_HALFWIDTH}\n"
                 f"(sirve un umbral 1D como el de ico?)")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.15)

    fig.suptitle(f"Tarea F.4 -- mapa (q4, q6) a nnn={nnn}, ideal {tname} "
                 f"(2e10 K/s, 300 K)")
    fig.tight_layout()
    out = f"taskF_mapa_nnn{nnn}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\n-> {out}")
