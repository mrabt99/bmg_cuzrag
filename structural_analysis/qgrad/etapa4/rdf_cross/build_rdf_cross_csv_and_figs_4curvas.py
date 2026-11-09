#!/usr/bin/env python3
"""Etapa 4 / Tarea 5 -- RDF cruzado Zr-Cu (g_ZrCu = g_CuZr, sanity de simetria incluido).
Tercer intento tras Tarea 3 (centro=especie, negativo) y Tarea 4 (mismo-tipo: Cu-Cu positivo,
Zr-Zr negativo). Completa la descomposicion binaria del RDF total con el par cruzado.

Entrada : raw/<tasa>/NN.rdf.txt (fix ave/time de c_rdf[*], 3 pares: total / Zr-Cu / Cu-Zr)
          y raw/<tasa>/NN.log (linea QC_BOX)
Salida  : rdf_cross_all.csv  rate,T,r,g_total,Ncum_total,g_ZrCu,Ncum_ZrCu,g_CuZr,Ncum_CuZr
          ../fig_K_rdf_heatmap_ZrCu.png, ../fig_L_rdf300_cross_zoom.png
Sanity  : g_total vs ../rdf/rdf_all.csv; g_ZrCu vs g_CuZr (deben coincidir, g_ij(r)=g_ji(r)).
T = 1900 - 100*(NN-1), NN=17 -> vidrio 300 K.
Sin pandas: numpy + csv.
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
RATES = ["1e12", "1e11", "2e10", "qgrad"]
RATE_LABEL = {"1e12": r"$10^{12}$ K/s", "1e11": r"$10^{11}$ K/s", "2e10": r"$2\times10^{10}$ K/s", "qgrad": "qgrad (graduada)"}
RATE_COLOR = {"1e12": "#d62728", "1e11": "#e6b800", "2e10": "#1f3fd6", "qgrad": "#7b3fbf"}
NBIN = 500


def T_of(nn):
    return 300 if nn == 17 else 1900 - 100 * (nn - 1)


def read_rdf(path):
    """7 valores por fila: r g_total coord_total g_ZrCu coord_ZrCu g_CuZr coord_CuZr."""
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) == 2:
                continue
            data.append([float(x) for x in p[1:8]])
    a = np.array(data)
    assert a.shape == (NBIN, 7), f"{path}: shape {a.shape}"
    return a


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
        r, g_tot, nc_tot, g_ZrCu, nc_ZrCu, g_CuZr, nc_CuZr = a.T
        if r_ref is None:
            r_ref = r
        assert np.allclose(r, r_ref), f"r-grid distinto en {rate} {nn}"
        T = T_of(nn)
        L, rho, nZr, nCu, nAg = read_qc(HERE / "raw" / rate / f"{nn:02d}.log")
        qc.append([rate, T, f"{nn:02d}", L, rho, nZr, nCu, nAg])
        for i in range(NBIN):
            tidy.append([rate, T, r[i], g_tot[i], nc_tot[i], g_ZrCu[i], nc_ZrCu[i], g_CuZr[i], nc_CuZr[i]])
        data[(rate, T)] = (r, g_tot, nc_tot, g_ZrCu, nc_ZrCu, g_CuZr, nc_CuZr)

with open(HERE / "rdf_cross_all.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rate", "T", "r", "g_total", "Ncum_total", "g_ZrCu", "Ncum_ZrCu", "g_CuZr", "Ncum_CuZr"])
    for row in tidy:
        w.writerow([row[0], row[1], f"{row[2]:.4f}"] + [f"{v:.6f}" for v in row[3:]])
print(f"rdf_cross_all.csv: {len(tidy)} filas (esperado {3*17*NBIN})")

comp = {(q[0], q[1]): (q[5], q[6], q[7]) for q in qc}
nZr0, nCu0, nAg0 = comp[(RATES[0], 300)]
for rate in RATES:
    for T in [300, 1900]:
        assert comp[(rate, T)] == (nZr0, nCu0, nAg0), f"composicion distinta en {rate} {T}"
print(f"composicion (fija): nZr={nZr0} nCu={nCu0} nAg={nAg0} (total {nZr0+nCu0+nAg0})")

# sanity 1: g_total vs Etapa 4 / Tarea 1
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
assert max_diff < 1e-6, "g_total del RDF cruzado no coincide con Etapa 4 / Tarea 1"

# sanity 2: g_ZrCu(r) == g_CuZr(r) (simetria del par cruzado, g_ij=g_ji)
max_asym = max(np.max(np.abs(data[(rate, T)][3] - data[(rate, T)][5]))
                for rate in RATES for T in {T_of(n) for n in range(1, 18)})
print(f"sanity simetria g_ZrCu vs g_CuZr: max|diff| = {max_asym:.2e}")
assert max_asym < 1e-6, "g_ZrCu(r) != g_CuZr(r), algo anda mal"

Ts = sorted({T_of(n) for n in range(1, 18)})
edges_T = np.concatenate([[Ts[0] - 50], (np.array(Ts[:-1]) + np.array(Ts[1:])) / 2, [Ts[-1] + 50]])
dr = r_ref[1] - r_ref[0]
edges_r = np.concatenate([r_ref - dr / 2, [r_ref[-1] + dr / 2]])
XMAX = 8.5

fig, axes = plt.subplots(4, 1, figsize=(6.2, 11.6), sharex=True)
norm = Normalize(vmin=0.0, vmax=3.0)
for ax, rate in zip(axes, RATES):
    Z = np.array([data[(rate, T)][3] for T in Ts])
    pm = ax.pcolormesh(edges_r, edges_T, Z, cmap="jet", norm=norm, shading="flat", rasterized=True)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(edges_T[0], edges_T[-1])
    ax.set_ylabel(r"$T$, K")
    ax.text(0.03, 0.93, RATE_LABEL[rate], transform=ax.transAxes, va="top", ha="left",
            fontsize=11, bbox=dict(fc="white", ec="0.3", lw=0.8, pad=2.5))
axes[-1].set_xlabel(r"$r$, $\mathrm{\AA}$")
axes[0].set_title("Par cruzado Zr-Cu", fontsize=11)
cax = fig.add_axes([0.87, 0.30, 0.025, 0.40])
cb = fig.colorbar(pm, cax=cax)
cb.set_label(r"$g_{\mathrm{Zr-Cu}}(r)$")
fig.subplots_adjust(left=0.13, right=0.84, top=0.96, bottom=0.07, hspace=0.06)
fig.savefig(OUT / "fig_K_rdf_heatmap_ZrCu_4curvas.png", dpi=200)
plt.close(fig)
print("fig_K_rdf_heatmap_ZrCu.png escrita")

fig, ax = plt.subplots(figsize=(5.2, 3.6))
for rate in RATES:
    r, g_tot, _, g_ZrCu, _, g_CuZr, _ = data[(rate, 300)]
    ax.plot(r, g_ZrCu, color=RATE_COLOR[rate], lw=1.6, label=RATE_LABEL[rate])
ax.set_xlim(3.0, 6.5)
ax.set_xlabel(r"$r$, $\mathrm{\AA}$")
ax.set_ylabel(r"$g_{\mathrm{Zr-Cu}}(r)$")
ax.set_title("glass 300 K -- par cruzado Zr-Cu", fontsize=10)
ax.legend(loc="upper right", frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "fig_L_rdf300_cross_zoom_4curvas.png", dpi=200)
plt.close(fig)
print("fig_L_rdf300_cross_zoom.png escrita")

# ---------------- Fig M: emergencia del hombro Zr-Cu vs T (analoga a fig_J de Cu-Cu) ----------------
from matplotlib.cm import ScalarMappable
Ts_show = [1900, 1500, 1100, 900, 700, 500, 300]
fig, ax = plt.subplots(figsize=(6.0, 5.0))
norm_T = Normalize(vmin=300, vmax=1900)
cmap_T = plt.get_cmap("viridis_r")
for T in Ts_show:
    r, _, _, g_ZrCu, _, _, _ = data[("qgrad", T)]
    ax.plot(r, g_ZrCu, color=cmap_T(norm_T(T)), lw=1.8, label=f"{T} K")
ax.set_xlim(3.5, 6.5)
ax.set_ylim(0.8, 2.2)
ax.set_xlabel(r"$r$, $\mathrm{\AA}$")
ax.set_ylabel(r"$g_{\mathrm{Zr-Cu}}(r)$")
ax.set_title(r"Par cruzado Zr-Cu vs T, perfil qgrad (1500-700 K graduada)")
sm = ScalarMappable(norm=norm_T, cmap=cmap_T)
cb = fig.colorbar(sm, ax=ax)
cb.set_label(r"$T$, K")
fig.tight_layout()
fig.savefig(OUT / "fig_M_ZrCu_cross_emergence_4curvas.png", dpi=200)
plt.close(fig)
print("fig_M_ZrCu_cross_emergence.png escrita")

# profundidad del valle del 2do pico. Ventanas AJUSTADAS a la forma real de g_ZrCu(r)
# (pico1 ~4.8-5.1, dip ~5.3-5.5, pico2 ~5.6-5.8 a 300K; a diferencia de Cu-Cu, que tiene
# pico1~4.4-4.5/dip~5.0/pico2~5.4 -- las ventanas de CuCu_split_summary.csv NO sirven aca).
rows = []
for rt in RATES:
    for T in sorted({T_of(n) for n in range(1, 18)}, reverse=True):
        r, g_tot, _, g_ZrCu, _, g_CuZr, _ = data[(rt, T)]
        m1 = (r > 4.4) & (r < 5.15); md = (r > 5.15) & (r < 5.55); m2 = (r > 5.5) & (r < 6.15)
        i1 = np.argmax(g_ZrCu[m1]); r1, g1 = r[m1][i1], g_ZrCu[m1][i1]
        idip = np.argmin(g_ZrCu[md]); rdip, gdip = r[md][idip], g_ZrCu[md][idip]
        i2 = np.argmax(g_ZrCu[m2]); r2v, g2 = r[m2][i2], g_ZrCu[m2][i2]
        depth = (g1 + g2) / 2 - gdip
        rows.append([rt, T, r1, g1, rdip, gdip, r2v, g2, depth])
with open(HERE / "ZrCu_split_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rate", "T", "r_peak1", "g_peak1", "r_dip", "g_dip", "r_peak2", "g_peak2", "depth"])
    for row in rows:
        w.writerow([row[0], row[1]] + [f"{v:.4f}" for v in row[2:]])
print(f"ZrCu_split_summary.csv: {len(rows)} filas")
