#!/usr/bin/env python3
"""Tarea F, punto 2 -- de donde sale el corte de q4 para Z16.

El criterio del paper (q4 <= 0.01 y |q6 - 0.1924| <= 0.01, nnn=16) da
CERO atomos en nuestras configs: el minimo real de q4(nnn=16) es 0.0129,
arriba del ideal tabulado 0.0108. La mitad de q6 si la pasan ~2383
atomos, asi que la estrategia es: etapa 1 = ventana de q6 (criterio del
paper tal cual), etapa 2 = corte de q4 SACADO DE LOS DATOS.

Este script NO fija el corte: muestra el histograma, busca si hay un gap
real y tabula que pasaria con cada candidato, para decidir mirando.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from boo_types import (SPECIES, SPECIES_COLOR, IDEAL, load_boo_dump,
                       features, mask_z16_q6window)

DUMP = "../boo/2e10/17.boo.dump"

data = load_boo_dump(DUMP)
types = data[:, 1].astype(int)
q4, q6, _w4, _w6 = features(data, 16)
ref = IDEAL["Z16"]

win = mask_z16_q6window(data)
q4w = q4[win]
tw = types[win]

print(f"=== {DUMP} ===")
print(f"Etapa 1 (|q6 - {ref['q6']}| <= 0.01, nnn=16): {win.sum()} atomos "
      f"({100.0*win.mean():.3f} % del sistema)")
for t, name in SPECIES.items():
    m = tw == t
    print(f"    {name:2s}: {m.sum():5d}  ({100.0*m.sum()/win.sum():5.2f} % del subconjunto, "
          f"{100.0*m.sum()/(types==t).sum():5.3f} % de esa especie)")

print(f"\n--- q4(nnn=16) del subconjunto vs de TODA la poblacion ---")
print(f"{'':12s} {'min':>8s} {'p1':>8s} {'p5':>8s} {'p10':>8s} {'p25':>8s} "
      f"{'p50':>8s} {'p90':>8s} {'max':>8s}")
for label, arr in (("subconjunto", q4w), ("todos", q4)):
    print(f"{label:12s} " + " ".join(
        f"{v:8.4f}" for v in (arr.min(), np.percentile(arr, 1),
                              np.percentile(arr, 5), np.percentile(arr, 10),
                              np.percentile(arr, 25), np.percentile(arr, 50),
                              np.percentile(arr, 90), arr.max())))
print(f"\nideal Z16 tabulado: q4 = {ref['q4']}  "
      f"(inalcanzable: minimo real del subconjunto = {q4w.min():.4f})")

# --- hay gap? mirar los saltos mas grandes en la cola de q4 bajo -------
print("\n--- busqueda de gap real en la cola de q4 baja del subconjunto ---")
s = np.sort(q4w)
ncheck = min(400, len(s) - 1)
gaps = s[1:ncheck + 1] - s[:ncheck]
order = np.argsort(gaps)[::-1][:8]
print(f"(8 saltos mas grandes entre atomos consecutivos, dentro de los "
      f"{ncheck} de q4 mas bajo)")
print(f"{'rank':>4s} {'salto':>9s} {'entre q4':>9s} {'y q4':>9s} "
      f"{'atomos por debajo':>18s}")
for r, i in enumerate(order, 1):
    print(f"{r:4d} {gaps[i]:9.5f} {s[i]:9.4f} {s[i+1]:9.4f} {i+1:18d}")
typical = np.median(np.diff(s[:ncheck + 1]))
print(f"salto tipico (mediana) en esa cola: {typical:.5f}  -> el mayor salto "
      f"es {gaps[order[0]]/typical:.1f}x el tipico")

# --- candidatos de corte ----------------------------------------------
print("\n--- que pasaria con cada candidato de corte de q4 ---")
cands = [
    ("p10 subconj.", np.percentile(q4w, 10)),
    ("p25 subconj.", np.percentile(q4w, 25)),
    ("p50 subconj.", np.percentile(q4w, 50)),
    ("p75 subconj.", np.percentile(q4w, 75)),
    ("p90 subconj.", np.percentile(q4w, 90)),
    ("p5 global",    np.percentile(q4,  5)),
    ("p10 global",   np.percentile(q4, 10)),
]
print(f"{'candidato':14s} {'q4_cut':>8s} {'N Z16':>7s} {'%sist':>7s} "
      + " ".join(f"{n:>7s}" for n in SPECIES.values()))
for label, cut in cands:
    sel = win & (q4 <= cut)
    row = f"{label:14s} {cut:8.4f} {sel.sum():7d} {100.0*sel.mean():7.3f} "
    row += " ".join(f"{(sel & (types == t)).sum():7d}" for t in SPECIES)
    print(row)

# --- histograma --------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))

ax = axes[0]
bins = np.linspace(0, max(q4w.max(), 0.16), 90)
ax.hist(q4, bins=bins, density=True, histtype="step", lw=1.4, color="#898781",
        label=f"todos los atomos (n={len(q4)})")
ax.hist(q4w, bins=bins, density=True, histtype="stepfilled", lw=1.6,
        color="#4a3aa7", alpha=0.25, label=f"ventana de q6 (n={win.sum()})")
ax.hist(q4w, bins=bins, density=True, histtype="step", lw=1.8, color="#4a3aa7")
ax.axvline(ref["q4"], color="#e34948", ls="--", lw=1.4)
ax.text(ref["q4"], ax.get_ylim()[1] * 0.97, f"  ideal Z16\n  q4={ref['q4']}",
        color="#e34948", fontsize=8, va="top")
ax.set_xlabel("q4 (nnn=16)")
ax.set_ylabel("densidad")
ax.set_title("q4 del subconjunto de la ventana de q6, vs toda la poblacion")
ax.legend(frameon=False, fontsize=8)
ax.grid(alpha=0.15)

ax = axes[1]
bins2 = np.linspace(0, 0.09, 70)
for t, name in SPECIES.items():
    m = tw == t
    if m.sum() == 0:
        continue
    ax.hist(q4w[m], bins=bins2, histtype="step", lw=1.8,
            color=SPECIES_COLOR[name], label=f"{name} (n={m.sum()})")
ax.axvline(ref["q4"], color="#e34948", ls="--", lw=1.4, label="ideal Z16")
for label, cut in (("p10", np.percentile(q4w, 10)),
                   ("p25", np.percentile(q4w, 25)),
                   ("p90", np.percentile(q4w, 90))):
    ax.axvline(cut, color="#0b0b0b", ls=":", lw=1.1)
    ax.text(cut, ax.get_ylim()[1] * 0.92, f" {label}", fontsize=7.5, rotation=90)
ax.set_xlabel("q4 (nnn=16)")
ax.set_ylabel("cuentas")
ax.set_title("zoom: q4 del subconjunto, por especie central")
ax.legend(frameon=False, fontsize=8)
ax.grid(alpha=0.15)

fig.suptitle("Tarea F.2 -- de donde sale el corte de q4 para Z16 "
             "(2e10 K/s, 300 K)")
fig.tight_layout()
fig.savefig("taskF_z16_q4.png", dpi=140)
plt.close(fig)
print("\n-> taskF_z16_q4.png")
