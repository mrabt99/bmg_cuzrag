#!/usr/bin/env python3
"""Localizacion cuantitativa de la deformacion en la traccion (OVITO). Uso:
    module load anaconda/2024.02-py311
    python localization.py --dir ../pilot_tension [--rate 2e10] [--thr 0.2] [--cutoff 3.75]
Lee <dir>/*.dump (cuadros cada d_eps = 0.005; cuadro k <-> eps = 0.005 k) y <dir>/*.curve.dat. Escribe en <dir>/analysis/localization/.
(1) eta de von Mises INCREMENTAL (AtomicStrainModifier, 'Shear Strain') en 10 ventanas consecutivas de d_eps = 0.02 entre eps = 0 y 0.20,
    con referencia en el cuadro de INICIO de cada ventana (no en eps = 0). affine_mapping = ToReference, eliminate_cell_deformation = False
    (incluye la deformacion homogenea de la ventana, ~0.03). Cutoff = primer minimo de la g(r) 300 K (rdf_all.csv, suavizado 5 bins) salvo --cutoff.
(2) psi = desvio estandar de eta entre atomos, por ventana (+ media, p99, fraccion con eta > thr). Figura: sigma-eps (arriba) y psi (abajo), eje eps COMPARTIDO.
(3) Render de la caja ENTERA solo con los atomos eta > thr, en las ventanas 0.06-0.10 y 0.16-0.20 (de 0.04, referencia en 0.06 y 0.16),
    vistas a lo largo de y (plano x-z) y de z (plano x-y).
(4) .npz con id y eta por ventana. Los IDs son los del restart del vidrio (se cruzan con la clasificacion estructural)."""
import sys, os, csv, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
from ovito.io import import_file
from ovito.modifiers import AtomicStrainModifier
from ovito.pipeline import ReferenceConfigurationModifier as R
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
D = os.path.abspath(arg("--dir", "../pilot_tension")); THR = float(arg("--thr", 0.2)); RATE = arg("--rate", "2e10")
RDF = "/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis/etapa4/rdf/rdf_all.csv"; DE = 0.005
dump = glob.glob(f"{D}/*.dump")[0]; curve = glob.glob(f"{D}/*.curve.dat")[0]; OUT = f"{D}/analysis/localization"; os.makedirs(OUT, exist_ok=True)
def first_min(rate, T="300", w=5):
    a = np.array([(float(r["r"]), float(r["g_total"])) for r in csv.DictReader(open(RDF)) if r["rate"] == rate and r["T"] in (T, T + ".0")])
    r, g = a[:, 0], np.convolve(a[:, 1], np.ones(w) / w, mode="same")
    ip = int(np.argmax(np.where((r > 2.3) & (r < 3.2), g, -1))); return r[int(np.argmin(np.where((r > r[ip]) & (r < 4.8), g, 1e9)))]
cutoff = float(arg("--cutoff", first_min(RATE))); print(f"{dump}\ncutoff = {cutoff:.2f} A ; umbral de render eta > {THR}")
def eta_inc(f0, f1):
    p = import_file(dump); p.modifiers.append(AtomicStrainModifier(cutoff=cutoff, reference_frame=f0, affine_mapping=R.AffineMapping.ToReference, eliminate_cell_deformation=False))
    d = p.compute(f1); o = np.argsort(d.particles["Particle Identifier"][...])
    return d, d.particles["Shear Strain"][...][o], d.particles["Particle Identifier"][...][o], d.particles.positions[...][o]
nfr = import_file(dump).source.num_frames; print("cuadros:", nfr)
# ---- (1)(2) 10 ventanas de 0.02 (4 cuadros)
W = [(4 * j, 4 * j + 4) for j in range(10)]; rows, ETA, IDS = [], [], None
for f0, f1 in W:
    d, eta, ids, _ = eta_inc(f0, f1); IDS = ids; ETA.append(eta.astype(np.float32))
    rows.append((f0 * DE, f1 * DE, 0.5 * (f0 + f1) * DE, eta.std(), eta.mean(), np.median(eta), np.percentile(eta, 99), (eta > THR).mean()))
    print(f"ventana {f0*DE:.3f}-{f1*DE:.3f}: psi = {eta.std():.4f}  media = {eta.mean():.4f}  mediana = {np.median(eta):.4f}  p99 = {np.percentile(eta,99):.3f}  frac(eta>{THR}) = {(eta>THR).mean():.4f}")
np.savez_compressed(f"{OUT}/eta_incremental_ventanas0.02.npz", id=IDS, eta=np.array(ETA), eps_ini=np.array([r[0] for r in rows]), eps_fin=np.array([r[1] for r in rows]), cutoff=cutoff)
with open(f"{OUT}/psi_ventanas.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["eps_ini", "eps_fin", "eps_centro", "psi_std_eta", "media_eta", "mediana_eta", "p99_eta", f"frac_eta>{THR}"]); w.writerows([[f"{x:.4f}" if i < 3 else f"{x:.5f}" for i, x in enumerate(r)] for r in rows])
c = np.loadtxt(curve, comments="#"); e, s = c[:, 1], c[:, 2]; ip = int(np.argmax(s)); R_ = np.array(rows)
fig, ax = plt.subplots(2, 1, figsize=(7.4, 6.4), sharex=True, gridspec_kw=dict(height_ratios=[1.35, 1]))
ax[0].plot(e, s, color="#2a78d6", lw=1.6); ax[0].set(ylabel="σxx (GPa)", ylim=(0, s.max() * 1.15), title=f"Localización de la deformación — {os.path.basename(D)}")
ax[1].plot(R_[:, 2], R_[:, 3], "o-", color="#0b0b0b", ms=5, lw=1.4, label="ψ = desvío estándar de η (ventanas de Δε = 0.02)")
ax[1].axhline(R_[0, 3], color="#8a8985", ls="--", lw=1.0, label="ψ de la primera ventana (0–0.02, casi elástica)")
ax[1].set(xlabel="deformación ingeniería ε (centro de la ventana)", ylabel="ψ (η incremental)", xlim=(0, 0.205), ylim=(0, None))
for a in ax:
    a.axvline(e[ip], color="#8a8985", ls=":", lw=1.0); a.grid(color="#e6e5e1", lw=.7); a.set_axisbelow(True)
    for sp in ("top", "right"): a.spines[sp].set_visible(False)
ax[0].text(e[ip] + 0.002, s.max() * 1.05, f"σmax, ε = {e[ip]:.3f}", fontsize=8.5, color="#52514e"); ax[1].legend(frameon=False, fontsize=8.5, loc="upper left")
fig.tight_layout(); fig.savefig(f"{OUT}/psi_vs_eps.png", dpi=150); plt.close(fig)
# ---- (3) renders de la caja entera, solo eta > THR
RW = [(12, 20), (32, 40)]; fig, axs = plt.subplots(2, 2, figsize=(12.5, 9), constrained_layout=True); out = {}
for r, (f0, f1) in enumerate(RW):
    d, eta, ids, pos = eta_inc(f0, f1); out[f"{f0*DE:.2f}-{f1*DE:.2f}"] = (ids, eta.astype(np.float32)); cell = np.array(d.cell[...]); lo = cell[:, 3]; L = np.diag(cell[:, :3]); m = eta > THR
    print(f"render {f0*DE:.2f}-{f1*DE:.2f} (Δε = {(f1-f0)*DE:.2f}): {int(m.sum())} atomos con eta > {THR} ({100*m.mean():.2f} %), psi = {eta.std():.4f}")
    for k, (a_, b_, nm) in enumerate([(0, 2, "a lo largo de y (plano x-z)"), (0, 1, "a lo largo de z (plano x-y)")]):
        ax = axs[r, k]; o = np.argsort(eta[m]); sc = ax.scatter(pos[m, a_][o], pos[m, b_][o], c=eta[m][o], s=5, cmap="viridis", vmin=THR, vmax=1.0, alpha=.75, linewidths=0)
        ax.plot([lo[a_], lo[a_] + L[a_], lo[a_] + L[a_], lo[a_], lo[a_]], [lo[b_], lo[b_], lo[b_] + L[b_], lo[b_] + L[b_], lo[b_]], color="#52514e", lw=.9)
        ax.set_aspect("equal"); ax.set_title(f"ε {f0*DE:.2f}→{f1*DE:.2f}, vista {nm}: {int(m.sum())} átomos con η > {THR}", fontsize=9.5)
        ax.set_xlabel("x (Å)"); ax.set_ylabel("z (Å)" if b_ == 2 else "y (Å)")
fig.colorbar(sc, ax=axs, shrink=.8, label="η incremental (von Mises)"); fig.savefig(f"{OUT}/render_eta_gt{THR}.png", dpi=140)
np.savez_compressed(f"{OUT}/eta_incremental_render_windows.npz", **{f"ids_{k}": v[0] for k, v in out.items()}, **{f"eta_{k}": v[1] for k, v in out.items()})
print("->", OUT)
