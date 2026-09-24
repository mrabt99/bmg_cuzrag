#!/usr/bin/env python3
"""Compuesto de las 3 tasas en UNA imagen (misma camara, misma escala, una leyenda con conteos).

Requiere haber corrido antes:  python ovito_snapshot.py --bare
Uso (desde el directorio que contiene ovito_out/):  python compose_snapshots.py
     con corte:  python ovito_snapshot.py --bare --slab 15 ; python compose_snapshots.py --slab 15
     (usa ovito_out/<tasa>_slab_bare.png y counts_slab_<tasa>.csv; escribe comparison_slab.png)
     variante despejada: python ovito_snapshot.py --bare --resto-transp 1 --radius-scale 0.65 ;
                         python compose_snapshots.py --resto-transp 1 --radius-scale 0.65   -> comparison_rt1_rs0.65.png
     por especie central: python ovito_snapshot.py --bare --resto-transp 1 --radius-scale 0.8 --center Cu ;
                          python compose_snapshots.py --resto-transp 1 --radius-scale 0.8 --center Cu   (idem Zr)
     sin ico-Voronoi:  python ovito_snapshot.py --bare --slab 15 --no-ico-voronoi ;
                       python compose_snapshots.py --slab 15 --no-ico-voronoi   -> comparison_slab_noicovor.png
Entrada : ovito_out/<tasa>_bare.png, ovito_out/counts_<tasa>.csv
Salida  : ovito_out/comparison_snapshots.png
Solo usa numpy + matplotlib (no importa ovito). Las 3 imagenes se recortan con el MISMO recuadro
(union de las zonas no blancas) para no alterar la escala relativa entre paneles.
"""
import csv
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = Path("ovito_out")
RATES = ["1e12", "1e11", "2e10"]              # rapido -> lento (igual orden que Fig A de la Etapa 4)
TITLE = {"2e10": r"$2\times10^{10}$ K/s", "1e11": r"$10^{11}$ K/s", "1e12": r"$10^{12}$ K/s"}
PAD = 12
args = sys.argv[1:]
SLAB = float(args[args.index("--slab") + 1]) if "--slab" in args else None
NOICO = "--no-ico-voronoi" in args
HIDE_RESTO = False
SUF = ("slab_" if SLAB else "") + ("noicovor_" if NOICO else "")
if "--resto-transp" in args or "--radius-scale" in args:         # mismo tag que ovito_snapshot.py
    rt = float(args[args.index("--resto-transp") + 1]) if "--resto-transp" in args else 0.5
    rs = float(args[args.index("--radius-scale") + 1]) if "--radius-scale" in args else 1.0
    SUF += f"rt{rt:g}_rs{rs:g}_"
    HIDE_RESTO = rt >= 1.0
CENTER = args[args.index("--center") + 1] if "--center" in args else None
if CENTER:
    SUF += f"{CENTER}_"
OUT_NAME = "comparison_" + (SUF.rstrip("_") if SUF else "snapshots") + ".png"

imgs = {r: plt.imread(OUT_DIR / f"{r}_{SUF}bare.png")[:, :, :3] for r in RATES}
assert len({im.shape for im in imgs.values()}) == 1, "las 3 imagenes deben tener el mismo tamano"
nonwhite = np.zeros(imgs[RATES[0]].shape[:2], bool)
for im in imgs.values():
    nonwhite |= (im < 0.97).any(axis=2)
rows, cols = np.where(nonwhite.any(axis=1))[0], np.where(nonwhite.any(axis=0))[0]
r0, r1 = max(rows[0] - PAD, 0), rows[-1] + PAD + 1
c0, c1 = max(cols[0] - PAD, 0), cols[-1] + PAD + 1

counts = {}
for r in RATES:
    with open(OUT_DIR / f"counts_{SUF}{r}.csv") as f:
        counts[r] = list(csv.DictReader(f))

h, w = r1 - r0, c1 - c0
fig = plt.figure(figsize=(15.0, 15.0 * (h / w) / 3 + 2.2))
gs = fig.add_gridspec(2, 3, height_ratios=[h / w * 1.0, 0.42], hspace=0.02, wspace=0.03,
                      left=0.01, right=0.99, top=(0.89 if CENTER else 0.95), bottom=0.01)
for j, r in enumerate(RATES):
    ax = fig.add_subplot(gs[0, j])
    ax.imshow(imgs[r][r0:r1, c0:c1])
    ax.set_axis_off()
    ax.set_title(TITLE[r] + ", 300 K" + (f" (slab {SLAB:g} Å)" if SLAB else "") + (f"\ncentros {CENTER} (% sobre átomos de {CENTER})" if CENTER else ""), fontsize=15)
    tx = fig.add_subplot(gs[1, j])
    tx.set_axis_off(); tx.set_xlim(0, 1); tx.set_ylim(0, 1)
    for i, row in enumerate(counts[r]):
        color = (float(row["r"]), float(row["g"]), float(row["b"]))
        n, ntot = int(row["N"]), int(row["n_total"])
        y = 0.90 - i * 0.16
        tx.scatter([0.20], [y], s=110, color=color, clip_on=False)
        nm = row["name"] + (" (no dibujado)" if row["name"] == "resto" and HIDE_RESTO else "")
        tx.text(0.25, y, nm, va="center", ha="left", fontsize=12, color=color)
        tx.text(0.97, y, f"{n}  ({100 * n / ntot:.1f}%)", va="center", ha="right", fontsize=12, color="0.15")
fig.savefig(OUT_DIR / OUT_NAME, dpi=150)
print("->", OUT_DIR / OUT_NAME, f"(recorte filas {r0}:{r1}, cols {c0}:{c1})")
