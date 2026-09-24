#!/usr/bin/env python3
"""Etapa 4 / Tarea 3 -- RDF parcial (centrado en Zr y en Cu), pedido de Jaroszewicz
tras ver la presentacion: probar si la separacion del 2do pico (Fig 5 Mendez, fingerprint
de vitrificacion) aparece al descomponer el g(r) total por especie central.

Entrada : raw/<tasa>/NN.rdf.txt (fix ave/time de c_rdf[*], 3 pares: total / Zr-centro / Cu-centro)
          y raw/<tasa>/NN.log (linea QC_BOX)
Salida  : rdf_partial_all.csv  rate,T,r,g_total,Ncum_total,g_Zr,Ncum_Zr,g_Cu,Ncum_Cu
          ../fig_D_rdf_heatmap_Zr.png, ../fig_E_rdf_heatmap_Cu.png, ../fig_F_rdf300_species_zoom.png
Sanity  : g_total debe coincidir con ../rdf/rdf_all.csv (misma cantidad, mismo r-grid, run 0
          sobre los mismos restarts).
T = 1900 - 100*(NN-1), NN=17 -> vidrio 300 K (igual que Etapa 2/3/4).
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
RATE_COLOR = {"1e12": "#d62728", "1e11": "#e6b800", "2e10": "#1f3fd6"}
NBIN = 500


def T_of(nn):
    return 300 if nn == 17 else 1900 - 100 * (nn - 1)


def read_rdf(path):
    """7 valores por fila: r g_total coord_total g_Zr coord_Zr g_Cu coord_Cu."""
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) == 2:            # "timestep nrows"
                continue
            data.append([float(x) for x in p[1:8]])
    a = np.array(data)
    assert a.shape == (NBIN, 7), f"{path}: shape {a.shape}"
    return a  # columns: r, g_tot, nc_tot, g_Zr, nc_Zr, g_Cu, nc_Cu


def read_qc(path):
    m = re.search(
        r"QC_BOX lx=(\S+) ly=(\S+) lz=(\S+) natoms=(\d+) numdens=(\S+) nZr=(\d+) nCu=(\d+) nAg=(\d+)",
        Path(path).read_text())
    lx, ly, lz, n, rho, nZr, nCu, nAg = m.groups()
    assert abs(float(lx) - float(ly)) < 1e-6 and abs(float(lx) - float(lz)) < 1e-6, "caja no cubica"
    assert int(n) == 65536
    return float(lx), float(rho), int(nZr), int(nCu), int(nAg)


tidy, qc, data = [], [], {}
r_ref = None
for rate in RATES:
    for nn in range(1, 18):
        a = read_rdf(HERE / "raw" / rate / f"{nn:02d}.rdf.txt")
        r, g_tot, nc_tot, g_Zr, nc_Zr, g_Cu, nc_Cu = a.T
        if r_ref is None:
            r_ref = r
        assert np.allclose(r, r_ref), f"r-grid distinto en {rate} {nn}"
        T = T_of(nn)
        L, rho, nZr, nCu, nAg = read_qc(HERE / "raw" / rate / f"{nn:02d}.log")
        qc.append([rate, T, f"{nn:02d}", L, rho, nZr, nCu, nAg])
        for i in range(NBIN):
            tidy.append([rate, T, r[i], g_tot[i], nc_tot[i], g_Zr[i], nc_Zr[i], g_Cu[i], nc_Cu[i]])
        data[(rate, T)] = (r, g_tot, nc_tot, g_Zr, nc_Zr, g_Cu, nc_Cu)

with open(HERE / "rdf_partial_all.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rate", "T", "r", "g_total", "Ncum_total", "g_Zr", "Ncum_Zr", "g_Cu", "Ncum_Cu"])
    for row in tidy:
        w.writerow([row[0], row[1], f"{row[2]:.4f}"] + [f"{v:.6f}" for v in row[3:]])
print(f"rdf_partial_all.csv: {len(tidy)} filas (esperado {3*17*NBIN})")

# composicion: sanity vs Etapa 4 (README: nCu=31130 mencionado para BADF)
comp = {(q[0], q[1]): (q[5], q[6], q[7]) for q in qc}
nZr0, nCu0, nAg0 = comp[(RATES[0], 300)]
for rate in RATES:
    for T in [300, 1900]:
        nZr, nCu, nAg = comp[(rate, T)]
        assert (nZr, nCu, nAg) == (nZr0, nCu0, nAg0), f"composicion distinta en {rate} {T}"
print(f"composicion (fija en las 51 configs): nZr={nZr0} nCu={nCu0} nAg={nAg0} (total {nZr0+nCu0+nAg0})")

# sanity: g_total vs Etapa 4 / Tarea 1 (rdf/rdf_all.csv), mismo run 0 sobre los mismos restarts
ref_path = HERE.parent / "rdf" / "rdf_all.csv"
ref = {}
with open(ref_path) as f:
    rd = csv.DictReader(f)
    for row in rd:
        ref.setdefault((row["rate"], row["T"]), []).append(float(row["g_total"]))
max_diff = 0.0
for rate in RATES:
    for nn in range(1, 18):
        T = T_of(nn)
        g_ref = np.array(ref[(rate, str(T))])
        g_new = data[(rate, T)][1]
        max_diff = max(max_diff, np.max(np.abs(g_ref - g_new)))
print(f"sanity g_total vs etapa4/rdf/rdf_all.csv: max|diff| = {max_diff:.2e}")
assert max_diff < 1e-6, "g_total del RDF parcial no coincide con Etapa 4 / Tarea 1"

# ---------------- Fig D / Fig E: heatmaps centrado-Zr / centrado-Cu ----------------
Ts = sorted({T_of(n) for n in range(1, 18)})
edges_T = np.concatenate([[Ts[0] - 50], (np.array(Ts[:-1]) + np.array(Ts[1:])) / 2, [Ts[-1] + 50]])
dr = r_ref[1] - r_ref[0]
edges_r = np.concatenate([r_ref - dr / 2, [r_ref[-1] + dr / 2]])
XMAX = 8.5

for species, idx, fname, title in [("Zr", 3, "fig_D_rdf_heatmap_Zr.png", r"$g_{\mathrm{Zr-all}}(r)$"),
                                    ("Cu", 5, "fig_E_rdf_heatmap_Cu.png", r"$g_{\mathrm{Cu-all}}(r)$")]:
    fig, axes = plt.subplots(3, 1, figsize=(6.2, 9.0), sharex=True)
    norm = Normalize(vmin=0.0, vmax=3.0)
    for ax, rate in zip(axes, RATES):
        Z = np.array([data[(rate, T)][idx] for T in Ts])
        pm = ax.pcolormesh(edges_r, edges_T, Z, cmap="jet", norm=norm, shading="flat", rasterized=True)
        ax.set_xlim(0, XMAX)
        ax.set_ylim(edges_T[0], edges_T[-1])
        ax.set_ylabel(r"$T$, K")
        ax.text(0.03, 0.93, RATE_LABEL[rate], transform=ax.transAxes, va="top", ha="left",
                fontsize=11, bbox=dict(fc="white", ec="0.3", lw=0.8, pad=2.5))
    axes[-1].set_xlabel(r"$r$, $\mathrm{\AA}$")
    axes[0].set_title(f"Centrado en {species}, vecino cualquier especie", fontsize=11)
    cax = fig.add_axes([0.87, 0.30, 0.025, 0.40])
    cb = fig.colorbar(pm, cax=cax)
    cb.set_label(title)
    fig.subplots_adjust(left=0.13, right=0.84, top=0.96, bottom=0.07, hspace=0.06)
    fig.savefig(OUT / fname, dpi=200)
    plt.close(fig)
    print(f"{fname} escrita")

# ---------------- Fig F: 300 K, total vs Zr-centrado vs Cu-centrado, zoom 2do pico ----------------
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), sharey=True)
for ax, idx, title in [(axes[0], 3, "centrado en Zr"), (axes[1], 5, "centrado en Cu")]:
    for rate in RATES:
        r, g_tot, _, g_Zr, _, g_Cu, _ = data[(rate, 300)]
        g = g_Zr if idx == 3 else g_Cu
        ax.plot(r, g, color=RATE_COLOR[rate], lw=1.6, label=RATE_LABEL[rate])
    ax.set_xlim(3.0, 6.5)
    ax.set_xlabel(r"$r$, $\mathrm{\AA}$")
    ax.set_title(f"glass 300 K -- {title}", fontsize=10)
    ax.axvspan(3.0, 3.0, alpha=0)  # no-op, keeps autoscale consistent
axes[0].set_ylabel(r"$g(r)$")
axes[0].legend(loc="upper right", frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "fig_F_rdf300_species_zoom.png", dpi=200)
plt.close(fig)
print("fig_F_rdf300_species_zoom.png escrita")
