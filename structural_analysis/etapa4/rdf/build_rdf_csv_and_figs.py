#!/usr/bin/env python3
"""Etapa 4 / Tarea 1 -- arma rdf_all.csv (tidy), rdf_qc.csv y Fig A / Fig B.

Entrada : raw/<tasa>/NN.rdf.txt (fix ave/time de c_rdf[*]) y raw/<tasa>/NN.log (linea QC_BOX)
Salida  : rdf_all.csv   rate,T,r,g_total,Ncum
          rdf_qc.csv    rate,T,config,L,numdens,Ncum_rmax,Ncum_ideal,ratio
          ../fig_A_rdf_heatmap.png, ../fig_B_rdf300_overlay.png
T = 1900 - 100*(NN-1), NN=17 -> glass 300 K (igual que Etapa 2/3).
Sin pandas (el venv no lo tiene): numpy + csv.
"""
import csv, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
RATES = ["1e12", "1e11", "2e10"]                    # orden de paneles: rapido arriba
RATE_LABEL = {"1e12": r"$10^{12}$ K/s", "1e11": r"$10^{11}$ K/s", "2e10": r"$2\times10^{10}$ K/s"}
RATE_COLOR = {"1e12": "#d62728", "1e11": "#e6b800", "2e10": "#1f3fd6"}  # rojo/amarillo/azul como Fig 3a
NBIN = 500


def T_of(nn):
    return 300 if nn == 17 else 1900 - 100 * (nn - 1)


def read_rdf(path):
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) == 2:            # "timestep nrows"
                continue
            data.append([float(x) for x in p[1:4]])
    a = np.array(data)
    assert a.shape == (NBIN, 3), f"{path}: shape {a.shape}"
    return a[:, 0], a[:, 1], a[:, 2]


def read_qc(path):
    m = re.search(r"QC_BOX lx=(\S+) ly=(\S+) lz=(\S+) natoms=(\d+) numdens=(\S+)", Path(path).read_text())
    lx, ly, lz, n, rho = m.groups()
    assert abs(float(lx) - float(ly)) < 1e-6 and abs(float(lx) - float(lz)) < 1e-6, "caja no cubica"
    assert int(n) == 65536
    return float(lx), float(rho)


tidy, qc, data = [], [], {}
r_ref = None
for rate in RATES:
    for nn in range(1, 18):
        r, g, nc = read_rdf(HERE / "raw" / rate / f"{nn:02d}.rdf.txt")
        if r_ref is None:
            r_ref = r
        assert np.allclose(r, r_ref), f"r-grid distinto en {rate} {nn}"
        T = T_of(nn)
        L, rho = read_qc(HERE / "raw" / rate / f"{nn:02d}.log")
        ideal = rho * 4 / 3 * np.pi * r[-1] ** 3       # N(<r) ideal en el ultimo bin (centro 9.99 A)
        qc.append([rate, T, f"{nn:02d}", L, rho, nc[-1], ideal, nc[-1] / ideal])
        for ri, gi, ci in zip(r, g, nc):
            tidy.append([rate, T, ri, gi, ci])
        data[(rate, T)] = (r, g, nc)

with open(HERE / "rdf_all.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rate", "T", "r", "g_total", "Ncum"])
    for row in tidy:
        w.writerow([row[0], row[1], f"{row[2]:.4f}", f"{row[3]:.6f}", f"{row[4]:.6f}"])
with open(HERE / "rdf_qc.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rate", "T", "config", "L", "numdens", "Ncum_rmax", "Ncum_ideal", "ratio"])
    for q in qc:
        w.writerow([q[0], q[1], q[2], f"{q[3]:.5f}", f"{q[4]:.6f}", f"{q[5]:.3f}", f"{q[6]:.3f}", f"{q[7]:.4f}"])

print(f"rdf_all.csv: {len(tidy)} filas (esperado {3*17*NBIN})")
ratios = np.array([q[7] for q in qc])
print(f"N(<rmax)/ideal: min={ratios.min():.4f} max={ratios.max():.4f}")
for rate in RATES:
    rho = [q[4] for q in qc if q[0] == rate]
    r300 = [q[4] for q in qc if q[0] == rate and q[1] == 300][0]
    print(f"  {rate}: numdens {min(rho):.5f}..{max(rho):.5f}  (300K: {r300:.5f})")

# ---------------- Fig A: heatmap g(r) en (r, T) ----------------
Ts = sorted({T_of(n) for n in range(1, 18)})
edges_T = np.concatenate([[Ts[0] - 50], (np.array(Ts[:-1]) + np.array(Ts[1:])) / 2, [Ts[-1] + 50]])
dr = r_ref[1] - r_ref[0]
edges_r = np.concatenate([r_ref - dr / 2, [r_ref[-1] + dr / 2]])
XMAX = 8.5
fig, axes = plt.subplots(3, 1, figsize=(6.2, 9.0), sharex=True)
norm = Normalize(vmin=0.0, vmax=3.0)
for ax, rate in zip(axes, RATES):
    Z = np.array([data[(rate, T)][1] for T in Ts])
    pm = ax.pcolormesh(edges_r, edges_T, Z, cmap="jet", norm=norm, shading="flat", rasterized=True)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(edges_T[0], edges_T[-1])
    ax.set_ylabel(r"$T$, K")
    ax.text(0.03, 0.93, RATE_LABEL[rate], transform=ax.transAxes, va="top", ha="left",
            fontsize=11, bbox=dict(fc="white", ec="0.3", lw=0.8, pad=2.5))
axes[-1].set_xlabel(r"$r$, $\mathrm{\AA}$")
cax = fig.add_axes([0.87, 0.30, 0.025, 0.40])
cb = fig.colorbar(pm, cax=cax)
cb.set_label(r"$g(r)$")
fig.subplots_adjust(left=0.13, right=0.84, top=0.985, bottom=0.07, hspace=0.06)
fig.savefig(OUT / "fig_A_rdf_heatmap.png", dpi=200)
plt.close(fig)

# ---------------- Fig B: g(r) a 300 K + N(<r) ----------------
fig, ax = plt.subplots(figsize=(5.2, 3.6))
ax2 = ax.twinx()
for rate in RATES:
    r, g, nc = data[(rate, 300)]
    ax.plot(r, g, color=RATE_COLOR[rate], lw=1.6, label=RATE_LABEL[rate])
    ax2.plot(r, nc, color=RATE_COLOR[rate], lw=1.3, ls=(0, (4, 2)))
ax.set_xlim(1.8, XMAX)
ax.set_ylim(0, 3.4)
ax.set_xlabel(r"$r$, $\mathrm{\AA}$")
ax.set_ylabel(r"$g(r)$")
ax2.set_ylim(0, 36)
ax2.set_yticks([12, 24])
ax2.set_ylabel(r"$N(<r)$")
ax.legend(loc="upper right", frameon=False, fontsize=9, title="glass 300 K", title_fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "fig_B_rdf300_overlay.png", dpi=200)
plt.close(fig)
print("figuras escritas")
