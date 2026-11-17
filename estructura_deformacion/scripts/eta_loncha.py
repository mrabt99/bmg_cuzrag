#!/usr/bin/env python3
"""Deformación de von Mises por átomo (OVITO AtomicStrain, η ACUMULADO desde ε = 0) para cruzarla con la estructura inicial
en los mapas de loncha de la bitácora (17/11/2026). Solo lee los dumps de tracción; escribe <caso>/analysis/eta_loncha.npz.

Parámetros (iguales a von_mises_strain.py del piloto): cutoff 3,75 Å (primer mínimo de g(r) a 300 K), referencia = cuadro 0
(ε = 0), affine_mapping = ToReference, eliminate_cell_deformation = False (η total).
Cuadros: ε = 0,02 / 0,04 / 0,06 / 0,08 / 0,10 / 0,15 / 0,20 (dump cada Δε = 0,005; el pico está en 0,055-0,068).
Los casos se buscan en mechanical/matrix/ (ruta absoluta), así que se puede correr desde cualquier lado.
Guarda: id, type, x0 (posiciones del cuadro 0, ordenadas por id), L0, eta_020 ... eta_200 (η acumulada hasta cada ε).
Uso (con OVITO): module load anaconda/2024.02-py311; python eta_loncha.py <caso> [<caso> ...]
"""
import sys, os, warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import numpy as np
from ovito.io import import_file
from ovito.modifiers import AtomicStrainModifier
from ovito.pipeline import ReferenceConfigurationModifier

CUTOFF = 3.75
FRAMES = {f"eta_{round(1000 * e):03d}": round(e / 0.005) for e in (0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20)}   # cuadro = ε/0,005
M = "/scratch/negishi/aaabt/bmg_cuzrag/mechanical/matrix"
for caso in (os.path.join(M, c) for c in sys.argv[1:]):
    dump = [f for f in os.listdir(caso) if f.endswith(".dump")][0]
    p = import_file(os.path.join(caso, dump))
    assert p.source.num_frames == 41, (caso, p.source.num_frames)
    d0 = p.compute(0)
    o = np.argsort(d0.particles["Particle Identifier"][...])
    out = dict(id=d0.particles["Particle Identifier"][...][o], type=d0.particles["Particle Type"][...][o],
               x0=np.array(d0.particles.positions[...])[o] - np.array(d0.cell[:, 3]), L0=np.diag(np.array(d0.cell[:, :3])))
    p.modifiers.append(AtomicStrainModifier(cutoff=CUTOFF, reference_frame=0, eliminate_cell_deformation=False,
                                            affine_mapping=ReferenceConfigurationModifier.AffineMapping.ToReference))
    for k, fr in FRAMES.items():
        d = p.compute(fr)
        oo = np.argsort(d.particles["Particle Identifier"][...])
        out[k] = d.particles["Shear Strain"][...][oo]
    np.savez_compressed(os.path.join(caso, "analysis", "eta_loncha.npz"), **out)
    print(os.path.basename(caso) + ": " + "  ".join(f"{k[4:]}: {out[k].mean():.3f}/{np.percentile(out[k], 99):.2f}" for k in FRAMES)
          + "   (η media/p99)", flush=True)
