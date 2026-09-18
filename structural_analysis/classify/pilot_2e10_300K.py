#!/usr/bin/env python3
"""Etapa 2b - PILOTO sobre UNA sola config antes de escalar a las 51.

Config elegida: 2e10 K/s, punto 17 (300 K relajado) -- la mas fria, donde
mas ordenamiento icosaedrico se espera.

Produce:
  1. distribucion de q6 (nnn=12) separada por especie central (Zr/Cu/Ag),
     para chequear el patron que describe Ryltsev: Cu-centrado icosaedrico,
     Zr-centrado mas Frank-Kasper/desordenado.
  2. histograma de distancias a cada tipo ideal (y de la distancia minima),
     por especie, para DECIDIR epsilon en vez de fijarlo a ojo.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from boo_types import (SPECIES, SPECIES_COLOR, IDEAL, DIST_TYPES,
                       EPSILON_DEFAULT, load_boo_dump, features,
                       distances_to_ideals, mask_ico, mask_z16)

DUMP = "../boo/2e10/17.boo.dump"

data = load_boo_dump(DUMP)
types = data[:, 1].astype(int)
natoms = data.shape[0]

print(f"=== {DUMP} : {natoms} atomos ===")
for t, name in SPECIES.items():
    n = (types == t).sum()
    print(f"  type {t} = {name:2s}: {n:6d}  ({100.0*n/natoms:5.2f} %)")

# ------------------------------------------------------------------ 1. q6
q4_12, q6_12, w4_12, w6_12 = features(data, 12)

print("\n=== q6 (nnn=12) por especie central ===")
print(f"{'especie':8s} {'media':>8s} {'mediana':>8s} {'p90':>8s} {'p99':>8s} "
      f"{'max':>8s} {'%q6>=0.6':>9s} {'%ico(q6&w6)':>12s}")
ico = mask_ico(data)
for t, name in SPECIES.items():
    m = types == t
    q6 = q6_12[m]
    pct_q6 = 100.0 * (q6 >= 0.6).sum() / m.sum()
    pct_ico = 100.0 * (ico & m).sum() / m.sum()
    print(f"{name:8s} {q6.mean():8.4f} {np.median(q6):8.4f} "
          f"{np.percentile(q6, 90):8.4f} {np.percentile(q6, 99):8.4f} "
          f"{q6.max():8.4f} {pct_q6:9.3f} {pct_ico:12.3f}")

print(f"\ntotal icosaedrico (criterio explicito del paper): "
      f"{ico.sum()} atomos = {100.0*ico.sum()/natoms:.3f} % del sistema")
z16 = mask_z16(data)
print(f"total Z16 (criterio explicito del paper): "
      f"{z16.sum()} atomos = {100.0*z16.sum()/natoms:.3f} % del sistema")
print("\nreparto de esos dos tipos por especie central:")
print(f"{'especie':8s} {'ico':>8s} {'% de esa especie':>18s} {'Z16':>8s} "
      f"{'% de esa especie':>18s}")
for t, name in SPECIES.items():
    m = types == t
    print(f"{name:8s} {(ico & m).sum():8d} {100.0*(ico & m).sum()/m.sum():18.3f} "
          f"{(z16 & m).sum():8d} {100.0*(z16 & m).sum()/m.sum():18.3f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
bins = np.linspace(0, 0.75, 120)
for t, name in SPECIES.items():
    m = types == t
    axes[0].hist(q6_12[m], bins=bins, density=True, histtype="step", lw=1.8,
                 color=SPECIES_COLOR[name], label=f"{name} (n={m.sum()})")
axes[0].axvline(0.6, color="#4a3aa7", ls="--", lw=1.2)
axes[0].text(0.6, axes[0].get_ylim()[1]*0.95, " q6=0.6\n (corte ico)",
             color="#4a3aa7", fontsize=8, va="top")
for t, name in [(1, "Zr"), (2, "Cu")]:
    axes[0].axvline(np.median(q6_12[types == t]), color=SPECIES_COLOR[name],
                    ls=":", lw=1.0)
axes[0].set_xlabel("q6 (nnn=12)")
axes[0].set_ylabel("densidad")
axes[0].set_title("2e10 K/s, 300 K: q6 por especie central")
axes[0].legend(frameon=False)
axes[0].grid(alpha=0.15)

# zoom a la cola icosaedrica
for t, name in SPECIES.items():
    m = types == t
    axes[1].hist(q6_12[m], bins=np.linspace(0.5, 0.72, 90), density=True,
                 histtype="step", lw=1.8, color=SPECIES_COLOR[name], label=name)
axes[1].axvline(0.6, color="#4a3aa7", ls="--", lw=1.2)
axes[1].set_xlabel("q6 (nnn=12)")
axes[1].set_ylabel("densidad")
axes[1].set_title("zoom a la cola (q6 > 0.5)")
axes[1].legend(frameon=False)
axes[1].grid(alpha=0.15)
fig.tight_layout()
fig.savefig("pilot_q6_por_especie.png", dpi=140)
plt.close(fig)
print("\n-> pilot_q6_por_especie.png")

# ------------------------------------------------- 2. distancias a ideales
dists = distances_to_ideals(data)          # (natoms, 4) para hcp/fcc/Z14/Z15
dmin = dists.min(axis=1)
argmin = dists.argmin(axis=1)

print("\n=== distancia a cada tipo ideal en (q4, q6, What6) ===")
print(f"{'tipo':6s} {'nnn':>4s} {'min':>8s} {'p1':>8s} {'p5':>8s} "
      f"{'mediana':>8s} {'%<0.05':>8s} {'%<0.03':>8s} {'%<0.10':>8s}")
for j, t in enumerate(DIST_TYPES):
    d = dists[:, j]
    print(f"{t:6s} {IDEAL[t]['nnn']:4d} {d.min():8.4f} "
          f"{np.percentile(d,1):8.4f} {np.percentile(d,5):8.4f} "
          f"{np.median(d):8.4f} {100.0*(d<0.05).mean():8.3f} "
          f"{100.0*(d<0.03).mean():8.3f} {100.0*(d<0.10).mean():8.3f}")

print("\n=== barrido de epsilon: % clasificado (solo tipos por distancia) ===")
print(f"{'epsilon':>8s} {'%asignado':>10s} " +
      " ".join(f"{t:>8s}" for t in DIST_TYPES) + f" {'%desord.':>9s}")
for eps in (0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.15):
    inside = dmin < eps
    row = f"{eps:8.2f} {100.0*inside.mean():10.3f} "
    for j, _t in enumerate(DIST_TYPES):
        row += f"{100.0*((argmin == j) & inside).mean():8.3f} "
    row += f"{100.0*(~inside).mean():9.3f}"
    print(row)

print(f"\n=== con epsilon = {EPSILON_DEFAULT}, reparto por especie central ===")
inside = dmin < EPSILON_DEFAULT
print(f"{'especie':8s} " + " ".join(f"{t:>8s}" for t in DIST_TYPES) +
      f" {'ico':>8s} {'Z16':>8s} {'desord.':>9s}")
for t, name in SPECIES.items():
    m = types == t
    row = f"{name:8s} "
    for j, _tt in enumerate(DIST_TYPES):
        row += f"{100.0*((argmin == j) & inside & m).sum()/m.sum():8.3f} "
    row += f"{100.0*(ico & m).sum()/m.sum():8.3f} "
    row += f"{100.0*(z16 & m).sum()/m.sum():8.3f} "
    row += f"{100.0*((~inside) & (~ico) & (~z16) & m).sum()/m.sum():9.3f}"
    print(row)

fig, axes = plt.subplots(2, 3, figsize=(15, 7.5))
axes = axes.ravel()
dbins = np.linspace(0, 0.5, 160)
for j, t in enumerate(DIST_TYPES):
    ax = axes[j]
    for tt, name in SPECIES.items():
        m = types == tt
        ax.hist(dists[m, j], bins=dbins, density=True, histtype="step", lw=1.6,
                color=SPECIES_COLOR[name], label=name)
    ax.axvline(EPSILON_DEFAULT, color="#4a3aa7", ls="--", lw=1.2)
    ax.set_title(f"distancia al ideal {t} (nnn={IDEAL[t]['nnn']})")
    ax.set_xlabel("distancia en (q4, q6, What6)")
    ax.set_ylabel("densidad")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.15)

ax = axes[4]
for tt, name in SPECIES.items():
    m = types == tt
    ax.hist(dmin[m], bins=dbins, density=True, histtype="step", lw=1.6,
            color=SPECIES_COLOR[name], label=name)
ax.axvline(EPSILON_DEFAULT, color="#4a3aa7", ls="--", lw=1.2)
ax.set_title("distancia MINIMA (al tipo mas cercano)")
ax.set_xlabel("distancia en (q4, q6, What6)")
ax.set_ylabel("densidad")
ax.legend(frameon=False, fontsize=8)
ax.grid(alpha=0.15)

ax = axes[5]
eps_grid = np.linspace(0.005, 0.2, 80)
for j, t in enumerate(DIST_TYPES):
    frac = [100.0 * ((argmin == j) & (dmin < e)).mean() for e in eps_grid]
    ax.plot(eps_grid, frac, lw=1.8, label=t)
ax.plot(eps_grid, [100.0 * (dmin < e).mean() for e in eps_grid], lw=2.0,
        color="#0b0b0b", ls="--", label="total asignado")
ax.axvline(EPSILON_DEFAULT, color="#4a3aa7", ls="--", lw=1.2)
ax.set_xlabel("epsilon")
ax.set_ylabel("% de atomos")
ax.set_title("poblacion asignada vs epsilon")
ax.legend(frameon=False, fontsize=8)
ax.grid(alpha=0.15)

fig.suptitle("2e10 K/s, 300 K -- distancias a los ideales de Ryltsev "
             "(linea punteada = epsilon 0.05 propuesto)")
fig.tight_layout()
fig.savefig("pilot_distancias.png", dpi=140)
plt.close(fig)
print("-> pilot_distancias.png")
