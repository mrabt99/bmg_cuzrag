#!/usr/bin/env python3
"""Tarea G, puntos 2a/2b/2c -- sobre el piloto ya calculado
(voronoi_pilot.py) para 2e10, config a 300K.

Convencion de columnas de 'Voronoi Index' en esta version de OVITO
(3.15.5), CONFIRMADA EMPIRICAMENTE (no solo por la doc, que resumida
por una herramienta de fetch resulto ambigua/inconsistente): la doc
dice "las 2 primeras componentes son siempre 0 porque un poligono
necesita >=3 aristas" e ilustra el icosaedro como (0,0,0,0,12,0,...).
Verificado a mano: 2992 atomos de esta config tienen el vector EXACTO
[0,0,0,0,12,0,0,0,0,0,0] con Coordination==12 (icosaedro Voronoi
estricto, un solo tipo de cara, las 12 pentagonales). Esa posicion es
la columna 4 (0-based). Por lo tanto columna i = caras con (i+1)
aristas: columna 2 = n3, columna 3 = n4, columna 4 = n5, columna 5 = n6.
El indice de Ryltsev/Voronoi <n3,n4,n5,n6> = columnas (2,3,4,5).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, "../classify")
from boo_types import load_boo_dump, mask_ico, SPECIES, SPECIES_COLOR  # noqa: E402

N3, N4, N5, N6 = 2, 3, 4, 5   # columnas de vindex para (n3,n4,n5,n6)

# ---------------------------------------------------------------- carga
vindex = np.load("vindex_2e10_300K.npy")
types_ovito = np.load("types_2e10_300K.npy")
coord = np.load("coord_2e10_300K.npy")

# necesitamos los Particle Identifier para cruzar contra BOO por id, no
# por orden de fila (no asumido que OVITO preserve el orden del dump)
from ovito.io import import_file
pipeline = import_file("../2e10/snaps/17.dump")
data = pipeline.compute(0)
ids_voro = data.particles["Particle Identifier"][...]
assert np.array_equal(types_ovito, data.particles["Particle Type"][...])

order = np.argsort(ids_voro)
ids_voro = ids_voro[order]
vindex = vindex[order]
types_voro = types_ovito[order]
coord = coord[order]
print(f"ids Voronoi: {ids_voro[0]}..{ids_voro[-1]}, ordenados y consecutivos: "
      f"{np.array_equal(ids_voro, np.arange(1, len(ids_voro)+1))}")

# ------------------------------------------------------------- punto 2a
strict_ico = (
    (vindex[:, N3] == 0) & (vindex[:, N4] == 0) &
    (vindex[:, N5] == 12) & (vindex[:, N6] == 0)
)
print(f"\n=== 2a: icosaedro Voronoi estricto <0,0,12,0> ===")
print(f"total: {strict_ico.sum()} atomos ({100.0*strict_ico.mean():.3f} % del sistema)")
print(f"{'especie':8s} {'N':>6s} {'% de esa especie':>18s}")
for t, name in SPECIES.items():
    m = types_voro == t
    n = (strict_ico & m).sum()
    print(f"{name:8s} {n:6d} {100.0*n/m.sum():18.3f}")

# ------------------------------------------------------------- punto 2b
boo_data = load_boo_dump("../boo/2e10/17.boo.dump")
ids_boo = boo_data[:, 0].astype(int)
order_boo = np.argsort(ids_boo)
ids_boo = ids_boo[order_boo]
boo_data = boo_data[order_boo]
types_boo = boo_data[:, 1].astype(int)

assert np.array_equal(ids_boo, ids_voro), (
    "los ids de boo/2e10/17.boo.dump y de 2e10/snaps/17.dump NO coinciden "
    "-- no son la misma config, no se puede cruzar")
assert np.array_equal(types_boo, types_voro), "tipos no coinciden entre boo y voro"

boo_ico = mask_ico(boo_data)   # criterio explicito del paper, nnn=12, ya validado

print(f"\n=== 2b: superposicion ico BOO (nnn=12, q6>=0.6 & What6<=-0.16) "
      f"vs ico Voronoi estricto <0,0,12,0> ===")
both = boo_ico & strict_ico
only_boo = boo_ico & ~strict_ico
only_voro = strict_ico & ~boo_ico
neither_but_flagged_total = boo_ico.sum() + strict_ico.sum() - both.sum()

print(f"{'':22s} {'total':>7s} " + " ".join(f"{n:>5s}" for n in SPECIES.values()))
def row(label, mask):
    counts = [int((mask & (types_voro == t)).sum()) for t in SPECIES]
    print(f"{label:22s} {mask.sum():7d} " + " ".join(f"{c:5d}" for c in counts))
row("BOO ico (total)", boo_ico)
row("Voronoi ico (total)", strict_ico)
row("coinciden (ambos)", both)
row("solo BOO", only_boo)
row("solo Voronoi", only_voro)
union = boo_ico | strict_ico
print(f"\nunion (marcado por al menos uno): {union.sum()} atomos")
print(f"jaccard (intersecc./union): {both.sum()}/{union.sum()} = "
      f"{both.sum()/union.sum():.3f}" if union.sum() else "sin datos")

# ------------------------------------------------------------- punto 2c
print(f"\n=== 2c: distribucion de numero de coordinacion por especie central ===")
print(f"{'especie':8s} {'media':>7s} {'mediana':>8s} {'moda':>6s} {'min':>5s} {'max':>5s}")
for t, name in SPECIES.items():
    m = types_voro == t
    cc = coord[m]
    vals, counts = np.unique(cc, return_counts=True)
    moda = vals[np.argmax(counts)]
    print(f"{name:8s} {cc.mean():7.3f} {np.median(cc):8.1f} {moda:6d} "
          f"{cc.min():5d} {cc.max():5d}")

fig, ax = plt.subplots(figsize=(7.5, 5))
cn_min, cn_max = coord.min(), coord.max()
bins = np.arange(cn_min - 0.5, cn_max + 1.5, 1)
width = 0.26
offsets = {1: -width, 2: 0, 3: width}
centers = np.arange(cn_min, cn_max + 1)
for t, name in SPECIES.items():
    m = types_voro == t
    cc = coord[m]
    counts = np.array([(cc == n).sum() for n in centers]) / m.sum() * 100.0
    ax.bar(centers + offsets[t], counts, width=width, color=SPECIES_COLOR[name],
           label=f"{name} (n={m.sum()})")
ax.set_xlabel("numero de coordinacion (Voronoi, CN)")
ax.set_ylabel("% de esa especie")
ax.set_title("2e10 K/s, 300 K -- CN por especie central (estilo Fig. 6, "
             "Mendez & Jaroszewicz)")
ax.set_xticks(centers)
ax.legend(frameon=False)
ax.grid(alpha=0.15, axis="y")
fig.tight_layout()
fig.savefig("cn_por_especie.png", dpi=140)
plt.close(fig)
print("\n-> cn_por_especie.png")
