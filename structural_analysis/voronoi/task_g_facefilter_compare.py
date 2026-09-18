#!/usr/bin/env python3
"""Tarea G, chequeo puntual -- repite 2a/2b con relative_face_threshold=0.01
y compara contra el resultado sin filtro (2992 / superposicion 97%,
jaccard 0.248, ver task_g_analysis.py)."""
import numpy as np
import sys

sys.path.insert(0, "../classify")
from boo_types import load_boo_dump, mask_ico, SPECIES  # noqa: E402

N3, N4, N5, N6 = 2, 3, 4, 5

def strict_ico_mask(vindex):
    return (
        (vindex[:, N3] == 0) & (vindex[:, N4] == 0) &
        (vindex[:, N5] == 12) & (vindex[:, N6] == 0)
    )

# --- confirmar que la convencion de columnas sigue siendo la misma con
# el filtro de area activo (no asumido, re-chequeado) -------------------
vindex_ff = np.load("vindex_2e10_300K_facefilter.npy")
coord_ff = np.load("coord_2e10_300K_facefilter.npy")
single12 = (vindex_ff == 12).any(axis=1) & ((vindex_ff > 0).sum(axis=1) == 1)
if single12.sum():
    col = np.argmax(vindex_ff[single12][0])
    print(f"chequeo columna: atomos con vector [.. ,12, ..] unico -> col {col} "
          f"(esperado 4, igual que sin filtro); coord de esos = "
          f"{coord_ff[single12][0].sum() if False else coord_ff[single12][:5]}")
    assert col == 4, "la convencion de columnas cambio con el filtro de area"

types_ff = np.load("types_2e10_300K_facefilter.npy")
ids_ff = np.load("ids_2e10_300K_facefilter.npy")
assert np.array_equal(ids_ff, np.arange(1, len(ids_ff) + 1)), "ids no consecutivos/ordenados"

# --- sin filtro (ya calculado antes) -----------------------------------
vindex_nf = np.load("vindex_2e10_300K.npy")
types_nf = np.load("types_2e10_300K.npy")
# reordenar por id tal como se hizo en task_g_analysis.py
from ovito.io import import_file
pipeline = import_file("../2e10/snaps/17.dump")
ids_nf = pipeline.compute(0).particles["Particle Identifier"][...]
order_nf = np.argsort(ids_nf)
vindex_nf = vindex_nf[order_nf]
types_nf = types_nf[order_nf]

strict_ico_nf = strict_ico_mask(vindex_nf)
strict_ico_ff = strict_ico_mask(vindex_ff)

print(f"\n=== 1. conteo de <0,0,12,0> con relative_face_threshold=0.01 ===")
print(f"{'especie':8s} {'sin filtro':>11s} {'con filtro 1%':>14s}")
for t, name in SPECIES.items():
    m = types_ff == t
    print(f"{name:8s} {(strict_ico_nf & (types_nf == t)).sum():11d} "
          f"{(strict_ico_ff & m).sum():14d}")
print(f"{'TOTAL':8s} {strict_ico_nf.sum():11d} {strict_ico_ff.sum():14d}")

# --- BOO, misma config -------------------------------------------------
boo_data = load_boo_dump("../boo/2e10/17.boo.dump")
ids_boo = boo_data[:, 0].astype(int)
order_boo = np.argsort(ids_boo)
boo_data = boo_data[order_boo]
ids_boo = ids_boo[order_boo]
assert np.array_equal(ids_boo, ids_ff), "ids de BOO no coinciden con los de voronoi/snaps"
boo_ico = mask_ico(boo_data)
types_boo = boo_data[:, 1].astype(int)
assert np.array_equal(types_boo, types_ff)

def overlap_table(voro_mask, label):
    both = boo_ico & voro_mask
    only_boo = boo_ico & ~voro_mask
    only_voro = voro_mask & ~boo_ico
    union = boo_ico | voro_mask
    jac = both.sum() / union.sum() if union.sum() else float("nan")
    print(f"\n=== 2. superposicion BOO vs Voronoi -- {label} ===")
    print(f"{'':22s} {'total':>7s} " + " ".join(f"{n:>5s}" for n in SPECIES.values()))
    def row(lbl, mask):
        counts = [int((mask & (types_ff == t)).sum()) for t in SPECIES]
        print(f"{lbl:22s} {mask.sum():7d} " + " ".join(f"{c:5d}" for c in counts))
    row("BOO ico (total)", boo_ico)
    row("Voronoi ico (total)", voro_mask)
    row("coinciden (ambos)", both)
    row("solo BOO", only_boo)
    row("solo Voronoi", only_voro)
    print(f"jaccard = {both.sum()}/{union.sum()} = {jac:.3f}")
    return jac, both.sum(), voro_mask.sum()

jac_ff, both_ff, total_ff = overlap_table(strict_ico_ff, "CON relative_face_threshold=0.01")

print(f"\n=== 3. comparacion filtro vs sin filtro ===")
print(f"{'':30s} {'sin filtro':>12s} {'con filtro 1%':>14s}")
print(f"{'N icosaedros Voronoi':30s} {strict_ico_nf.sum():12d} {strict_ico_ff.sum():14d}")
print(f"{'coinciden con BOO (748 antes)':30s} {748:12d} {both_ff:14d}")
print(f"{'jaccard':30s} {0.248:12.3f} {jac_ff:14.3f}")
print(f"{'razon Voronoi/BOO (770)':30s} {strict_ico_nf.sum()/770:12.2f} "
      f"{strict_ico_ff.sum()/770:14.2f}")
