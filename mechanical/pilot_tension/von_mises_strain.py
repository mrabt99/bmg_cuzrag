#!/usr/bin/env python3
"""Deformacion atomica de von Mises (OVITO, AtomicStrainModifier -> propiedad 'Shear Strain' = eta) del cuadro final (strain 0.20) respecto del
cuadro 0 (strain 0) del dump de la traccion. Necesita OVITO: usar el python de  module load anaconda/2024.02-py311  (el venv_structural no lo tiene).
Cutoff = primer minimo de la g(r) total del vidrio 2e10 a 300 K (structural_analysis/etapa4/rdf/rdf_all.csv, SOLO LECTURA): maximo del primer pico entre 2.3 y 3.2 A,
minimo siguiente hasta 4.8 A, sobre g suavizada con media movil de 5 bins (el g viene de UN snapshot y es ruidoso). Se puede forzar con --cutoff.
Opciones de OVITO: cuadro de referencia 0, affine_mapping = ToReference (la caja cambia de forma), eliminate_cell_deformation = False -> eta es la deformacion TOTAL
(su media debe parecerse al valor homogeneo de Green-Lagrange calculado con las dimensiones de la caja; se imprime como verificacion).
Salidas (analysis/): eta_hist_eps0.20.csv/.png (histograma; en gris el de strain 0 -> 0.005 como referencia de ruido termico), eta_slice_eps0.20.png (cortes finos), eta_resumen.txt, eta_eps0.20.npz (id, eta).
Uso: python von_mises_strain.py [dump] [--cutoff X] [--frame N] [--slab A]"""
import sys, os, csv, warnings
warnings.filterwarnings("ignore")
import numpy as np
from ovito.io import import_file
from ovito.modifiers import AtomicStrainModifier
from ovito.pipeline import ReferenceConfigurationModifier as R
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
dump = next((a for a in sys.argv[1:] if a.endswith(".dump")), "tens_2e10_300K_1e9_x.dump")
RDF = "/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis/etapa4/rdf/rdf_all.csv"
SLAB = float(arg("--slab", 10.0)); os.makedirs("analysis", exist_ok=True)

def first_min(rate="2e10", T="300", w=5):
    a = np.array([(float(r["r"]), float(r["g_total"])) for r in csv.DictReader(open(RDF)) if r["rate"] == rate and r["T"] in (T, T + ".0")])
    r, g = a[:, 0], np.convolve(a[:, 1], np.ones(w) / w, mode="same")
    ip = int(np.argmax(np.where((r > 2.3) & (r < 3.2), g, -1))); im = int(np.argmin(np.where((r > r[ip]) & (r < 4.8), g, 1e9)))
    return r[ip], r[im], g[im]
rp, rm, gm = first_min(); cutoff = float(arg("--cutoff", rm))
print(f"g(r) 2e10 300 K: primer pico r = {rp:.2f} A ; primer minimo r = {rm:.2f} A (g = {gm:.3f}) -> cutoff = {cutoff:.2f} A")

pipe = import_file(dump); nf = pipe.source.num_frames; fin = int(arg("--frame", nf - 1)); print(f"{dump}: {nf} cuadros; se calcula el cuadro {fin} respecto del 0")
def eta_at(frame, ref=0):
    p = import_file(dump); p.modifiers.append(AtomicStrainModifier(cutoff=cutoff, reference_frame=ref, affine_mapping=R.AffineMapping.ToReference, eliminate_cell_deformation=False))
    d = p.compute(frame); return d, d.particles["Shear Strain"][...], d.particles["Particle Identifier"][...]
d, eta, ids = eta_at(fin); d1, eta1, _ = eta_at(1); d0, eta0, _ = eta_at(0)
c0 = np.array(import_file(dump).compute(0).cell[...]); cN = np.array(d.cell[...]); L0, LN = np.diag(c0[:, :3]), np.diag(cN[:, :3])
F = np.diag(LN / L0); E = 0.5 * (F.T @ F - np.eye(3)); ex, ey, ez = np.diag(E); vm_h = np.sqrt(0.5 * ((ex - ey) ** 2 + (ey - ez) ** 2 + (ez - ex) ** 2))
eps_x = LN[0] / L0[0] - 1
q = lambda x, p: np.percentile(x, p)
lines = [f"cutoff = {cutoff:.2f} A ; cuadro {fin} (strain ingenieria en x = {eps_x:.4f}) respecto del cuadro 0 ; N = {len(eta)}",
         f"verificacion cuadro 0 vs 0: max eta = {eta0.max():.2e} (debe ser 0)",
         f"eta von Mises (cuadro {fin}): media = {eta.mean():.4f}  mediana = {np.median(eta):.4f}  desvio = {eta.std():.4f}  p5/p95/p99/max = {q(eta,5):.3f}/{q(eta,95):.3f}/{q(eta,99):.3f}/{eta.max():.3f}",
         f"valor HOMOGENEO esperado (von Mises de Green-Lagrange con L/L0 = {np.round(LN/L0,4)}) = {vm_h:.4f}   (la media de eta deberia rondarlo: {eta.mean():.4f})",
         f"fraccion de atomos con eta > 2x mediana = {(eta > 2*np.median(eta)).mean():.4f} ; eta > 0.5 = {(eta > 0.5).mean():.4f} ; eta > 1.0 = {(eta > 1.0).mean():.5f}",
         f"referencia de ruido: cuadro 1 (strain {import_file(dump).compute(1).cell[...][0,0]/L0[0]-1:.4f}) vs 0: media = {eta1.mean():.4f}  mediana = {np.median(eta1):.4f}  p99 = {q(eta1,99):.3f}"]
print("\n".join(lines)); open("analysis/eta_resumen.txt", "w").write("\n".join(lines) + "\n")
o = np.argsort(ids); np.savez_compressed("analysis/eta_eps0.20.npz", id=ids[o], eta=eta[o])

# ---- histograma
bins = np.linspace(0, max(1.8, eta.max()), 91); h, ed = np.histogram(eta, bins=bins, density=True); h1, _ = np.histogram(eta1, bins=bins, density=True)
with open("analysis/eta_hist_eps0.20.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["eta_lo", "eta_hi", "densidad_eps0.20", "densidad_eps0.005_referencia", "conteo_eps0.20"]); cnt = np.histogram(eta, bins=bins)[0]
    for k in range(len(h)): w.writerow([f"{ed[k]:.4f}", f"{ed[k+1]:.4f}", f"{h[k]:.5f}", f"{h1[k]:.5f}", cnt[k]])
fig, ax = plt.subplots(figsize=(6.6, 4.2)); c = 0.5 * (ed[1:] + ed[:-1])
ax.fill_between(c, h1, color="#c9c8c3", alpha=.8, label=f"ε 0 → 0.005 (referencia de ruido)"); ax.plot(c, h, color="#2a78d6", lw=1.8, label=f"ε 0 → {eps_x:.2f}")
ax.axvline(eta.mean(), color="#0b0b0b", lw=.9, ls=":"); ax.text(eta.mean() + .02, ax.get_ylim()[1] * .92, f"media {eta.mean():.2f}", fontsize=8.5)
ax.set(xlabel="η von Mises (cutoff %.2f Å)" % cutoff, ylabel="densidad de probabilidad", xlim=(0, bins[-1]), title="Deformación atómica de von Mises — 2e10, 300 K, 1e9 s⁻¹"); ax.legend(frameon=False, fontsize=8.5)
ax.grid(color="#e6e5e1", lw=.7); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("analysis/eta_hist_eps0.20.png", dpi=150); plt.close(fig)

# ---- cortes finos: caja del cuadro final, coloreado por eta (escala comun 0 - 1.0)
pos = d.particles.positions[...]; lo = cN[:, 3]; ex_ = np.diag(cN[:, :3]); VMAX = 1.0
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
for ax, (ax_i, a_, b_, tt) in zip(axs, [(2, 0, 1, "corte en z (vista x-y)"), (1, 0, 2, "corte en y (vista x-z)")]):
    mid = lo[ax_i] + ex_[ax_i] / 2; s = np.abs(pos[:, ax_i] - mid) < SLAB / 2; oo = np.argsort(eta[s])
    sc = ax.scatter(pos[s, a_][oo], pos[s, b_][oo], c=eta[s][oo], s=3.2, cmap="viridis", vmin=0, vmax=VMAX, linewidths=0)
    ax.set_aspect("equal"); ax.set_title(f"{tt}, espesor {SLAB:.0f} Å ({int(s.sum())} átomos)", fontsize=10); ax.set_xlabel("x (Å) — eje de tracción →"); ax.set_ylabel("y (Å)" if b_ == 1 else "z (Å)")
fig.colorbar(sc, ax=axs, shrink=.85, label=f"η von Mises (saturado en {VMAX:.1f})")
fig.suptitle(f"Deformación atómica de von Mises, ε = {eps_x:.2f} respecto de ε = 0 (2e10, 300 K, 1e9 s⁻¹)", fontsize=11); fig.savefig("analysis/eta_slice_eps0.20.png", dpi=150)
print("-> analysis/: eta_resumen.txt, eta_hist_eps0.20.{csv,png}, eta_slice_eps0.20.png, eta_eps0.20.npz")
