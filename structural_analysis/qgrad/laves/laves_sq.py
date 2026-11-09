#!/usr/bin/env python3
"""Firma de largo alcance (e): S(k) 3D por FFT (asignacion NGP a una grilla 192^3; el factor de ventana es igual en todas las configs,
solo se comparan entre si) sobre vidrios a 300 K (2e10, 1e11, 1e12, qgrad 17u) y sobre las esferas implantadas de C15 (calibracion).
Para cada |k| (bines de 0.02 A^-1, 0.6-3.6 A^-1) se guarda media y MAXIMO de S(k)=|rho_k|^2/N. Un grano cristalino da picos aislados de S(k)
(max >> media a los q de Bragg del C15: (111)=1.53, (220)=2.50, (311)=2.93, (222)=3.06 A^-1 para a=7.10 A); un vidrio da max/media ~ log(N_k).
Salida: laves_sq.csv (max/media por config en las bandas de Bragg), fig_laves_sq.png"""
import csv, glob, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SA = "/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis/"; NG = 192
def read_pos(path):
    lines = open(path).read().split("\n"); n = int(lines[3]); b = np.array([l.split()[:2] for l in lines[5:8]], float)
    a = np.array([l.split() for l in lines[9:9 + n]], float); return a[:, 2:5] - b[:, 0], (b[:, 1] - b[:, 0])[0]
def sk(path):
    x, L = read_pos(path); n = len(x)
    g = np.floor(x / L * NG).astype(int) % NG
    rho = np.bincount((g[:, 0] * NG + g[:, 1]) * NG + g[:, 2], minlength=NG ** 3).reshape(NG, NG, NG).astype(float)
    F = np.fft.rfftn(rho); S = np.abs(F) ** 2 / n
    kx = np.fft.fftfreq(NG, d=1.0 / NG) * 2 * np.pi / L; kz = np.fft.rfftfreq(NG, d=1.0 / NG) * 2 * np.pi / L
    K = np.sqrt(kx[:, None, None] ** 2 + kx[None, :, None] ** 2 + kz[None, None, :] ** 2)
    w = np.ones_like(S); w[:, :, 1:-1] = 2.0                        # rfft: mitad del espacio k
    return K.ravel(), S.ravel(), w.ravel()
bins = np.arange(0.6, 3.62, 0.02); qc = 0.5 * (bins[1:] + bins[:-1])
def profile(path):
    K, S, w = sk(path); i = np.digitize(K, bins) - 1; ok = (i >= 0) & (i < len(qc))
    mean = np.bincount(i[ok], weights=(S * w)[ok], minlength=len(qc))[:len(qc)] / np.bincount(i[ok], weights=w[ok], minlength=len(qc))[:len(qc)]
    mx = np.zeros(len(qc)); np.maximum.at(mx, i[ok], S[ok]); return mean, mx
BRAGG = {"(111)": 1.53, "(220)": 2.50, "(311)": 2.93, "(222)": 3.06}
files = {"2e10 300K": SA + "2e10/snaps/17.dump", "1e11 300K": SA + "1e11/snaps/17.dump", "1e12 300K": SA + "1e12/snaps/17.dump", "qgrad 300K (%s)" % ("17r" if os.path.exists(SA + "qgrad/snaps/17r.dump") else "17u"): SA + "qgrad/snaps/%s.dump" % ("17r" if os.path.exists(SA + "qgrad/snaps/17r.dump") else "17u")}
for R in (8, 10, 15, 20): files[f"implante C15 R={R} A"] = f"implant/C15_R{R}.dump"
res = {k: profile(p) for k, p in files.items()}
rows = []
print(f"{'config':24s}  max S en [Bragg +-0.06] (media local del vidrio en parentesis) ; y max global S(k) en 1.0-3.4")
for k, (mean, mx) in res.items():
    row = [k]
    for nm, q in BRAGG.items():
        m = np.abs(qc - q) <= 0.06; row += [round(float(mx[m].max()), 2)]
    g = (qc > 1.0) & (qc < 3.4); row += [round(float(mx[g].max()), 2)]
    rows.append(row); print(f"{k:24s} " + "  ".join(f"{nm}={v:7.2f}" for nm, v in zip(list(BRAGG) + ['max_global'], row[1:])))
with open("laves_sq.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["config"] + [f"maxS_{n}" for n in BRAGG] + ["maxS_global_1.0-3.4"]); w.writerows(rows)
fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
for k, (mean, mx) in res.items():
    st = "--" if "implante" in k else "-"
    ax[0].plot(qc, mx, st, lw=1, label=k)
for q in BRAGG.values(): ax[0].axvline(q, color="grey", lw=.6, ls=":")
ax[0].set(xlabel="k (1/A)", ylabel="max S(k) en cada bin de |k|", yscale="log", title="Maximo de S(k) por |k| (lineas: Bragg C15)"); ax[0].legend(fontsize=7)
for k in [k for k in res if k.startswith(("2e10", "qgrad"))]: ax[1].plot(qc, res[k][0], lw=1.2, label=k)
ax[1].set(xlabel="k (1/A)", ylabel="<S(k)> (promedio angular)", title="S(k) promedio: qgrad vs 2e10"); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig("fig_laves_sq.png", dpi=150)
