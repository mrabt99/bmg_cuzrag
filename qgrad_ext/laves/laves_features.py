#!/usr/bin/env python3
"""Descriptores por atomo para la busqueda de Laves (Cu2Zr C15/C14). Corre con anaconda (OVITO). Solo lectura de snaps/ y boo/.
Por configuracion guarda features/<set>/<label>.npz con (ordenado por id):
  type, vindex(n3..n6), coord (Voronoi DEFAULT de OVITO, igual que el resto del pipeline),
  nZr (n de vecinos Voronoi que son Zr), tdev (solo Zr con <0,0,12,4> y exactamente 4 vecinos Zr: max |angulo - 109.47| entre
  los 6 pares de esos 4 vecinos, grados; NaN si no aplica), diamond (Zr: 0 otro, 1 cubico, 4 hexagonal, ...; IdentifyDiamond
  sobre la SUBRED de Zr, Cu/Ag borrados; resto 0), q4,q6,w6 nnn=12 y q4,q6,w6 nnn=16 (BOO, dumps de validate_boo.lmp),
  neigh_ptr/neigh_idx (grafo de vecinos Voronoi, CSR, para los clusters).
Uso: python laves_features.py <set> <label> <pos_dump> <boo_dump> [frame]
  set: refs | 2e10 | 1e11 | 1e12 | qgrad ; frame: para dumps multi-frame (refs), default ultimo."""
import sys, warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import numpy as np
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier, SelectTypeModifier, DeleteSelectedModifier, IdentifyDiamondModifier
TET = 109.4712

def read_boo_frame(path, frame):
    frames, cur = [], None
    with open(path) as f:
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                if cur is not None: frames.append(cur)
                cur = []; continue
            if line.startswith("ITEM:"): continue
            if cur is not None and line.strip() and len(line.split()) > 10: cur.append(line.split())
    frames.append(cur)
    a = np.array(frames[frame], float)
    return a[np.argsort(a[:, 0])]

def compute(pos_dump, boo_dump, frame):
    p = import_file(pos_dump)
    if frame is None: frame = p.source.num_frames - 1
    q = import_file(pos_dump)
    p.modifiers.append(VoronoiAnalysisModifier(compute_indices=True, generate_bonds=True))     # DEFAULT
    d = p.compute(frame)
    ids = d.particles["Particle Identifier"][...]; o = np.argsort(ids)
    inv = np.empty_like(o); inv[o] = np.arange(len(o))
    t = d.particles["Particle Type"][...][o]; vi = d.particles["Voronoi Index"][...][o][:, 2:6]; cn = d.particles["Coordination"][...][o]
    x = d.particles.positions[...][o]; Hm = np.array(d.cell[...])[:, :3]
    b = d.particles.bonds; top = inv[b.topology[...]]; pi = b["Periodic Image"][...]
    vec = x[top[:, 1]] - x[top[:, 0]] + pi @ Hm.T                       # vector a->b
    n = len(t)
    src = np.concatenate([top[:, 0], top[:, 1]]); dst = np.concatenate([top[:, 1], top[:, 0]]); V = np.concatenate([vec, -vec])
    order = np.argsort(src, kind="stable"); src, dst, V = src[order], dst[order], V[order]
    ptr = np.concatenate([[0], np.cumsum(np.bincount(src, minlength=n))])
    nZr = np.bincount(src, weights=(t[dst] == 1), minlength=n).astype(int)
    tdev = np.full(n, np.nan)
    cand = np.where((t == 1) & (vi == np.array([0, 0, 12, 4])).all(1) & (nZr == 4))[0]
    iu = np.triu_indices(4, 1)
    for i in cand:
        sl = slice(ptr[i], ptr[i + 1]); m = t[dst[sl]] == 1; u = V[sl][m]; u = u / np.linalg.norm(u, axis=1, keepdims=True)
        ang = np.degrees(np.arccos(np.clip((u @ u.T)[iu], -1, 1))); tdev[i] = np.abs(ang - TET).max()
    # subred Zr -> diamante
    r = import_file(pos_dump)
    r.modifiers.append(SelectTypeModifier(property="Particle Type", types={2, 3})); r.modifiers.append(DeleteSelectedModifier())
    r.modifiers.append(IdentifyDiamondModifier()); dd = r.compute(frame)
    idz = dd.particles["Particle Identifier"][...]; st = dd.particles["Structure Type"][...]
    diamond = np.zeros(n, int); pos_of_id = {int(k): i for i, k in enumerate(ids[o])}
    for k, s in zip(idz, st): diamond[pos_of_id[int(k)]] = s
    if boo_dump is None:                      # pruebas sinteticas (implantes): sin BOO
        B = np.full((n, 18), np.nan)
    else:
        B = read_boo_frame(boo_dump, frame); assert len(B) == n and np.array_equal(B[:, 0].astype(int), ids[o]), "orden de ids BOO != posiciones"
        assert np.array_equal(B[:, 1].astype(int), t)
    # columnas del dump BOO: 2..5 = Q4 Q6 W4 W6 (nnn12), 14..17 = (nnn16)
    return dict(type=t, vindex=vi, coord=cn, nZr=nZr, tdev=tdev, diamond=diamond,
                q4_12=B[:, 2], q6_12=B[:, 3], w6_12=B[:, 5], q4_16=B[:, 14], q6_16=B[:, 15], w6_16=B[:, 17],
                neigh_ptr=ptr, neigh_idx=dst)

if __name__ == "__main__":
    st, lab, pos, boo = sys.argv[1:5]; fr = int(sys.argv[5]) if len(sys.argv) > 5 else None
    import os; os.makedirs(f"features/{st}", exist_ok=True)
    F = compute(pos, boo, fr); np.savez_compressed(f"features/{st}/{lab}.npz", **F)
    t = F["type"]; z16 = (t == 1) & (F["vindex"] == np.array([0, 0, 12, 4])).all(1)
    print(f"{st}/{lab}: Zr-Z16={int(z16.sum())}  Zr-Z16 con 4 Zr-vecinos={int(np.isfinite(F['tdev']).sum())}  diamante(Zr)={int((F['diamond']>0).sum())}", flush=True)
