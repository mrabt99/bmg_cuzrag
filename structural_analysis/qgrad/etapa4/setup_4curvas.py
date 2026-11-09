#!/usr/bin/env python3
"""Etapa 4 para qgrad: arma qgrad/etapa4/<tipo>/ con
  (a) raw/{2e10,1e11,1e12} = SYMLINKS a etapa4/<tipo>/raw/<tasa> (solo lectura, no se duplica ni se toca nada),
  (b) raw/qgrad/NN.{rdf.txt,log}: NN 01-05 = symlinks a los raw de 2e10 (mismo restart, no se recalcula),
      06-16 = ya calculados (etiquetas 06..16, que coinciden con NN), 17 = symlink a 17r (relajado) si existe, si no a 17u,
  (c) copia PARCHEADA (sufijo _4curvas) de cada builder original: RATES + qgrad, 4 paneles en los heatmaps, nombres de
      salida con sufijo _4curvas. Los builders/figuras/CSV originales no se tocan.
Uso: python setup_4curvas.py   (desde qgrad/etapa4/; idempotente)"""
import os, re
from pathlib import Path
HERE = Path(__file__).resolve().parent
E4 = HERE.parent.parent / "etapa4"                 # structural_analysis/etapa4 (originales, SOLO LECTURA)
KINDS = {"rdf": "build_rdf_csv_and_figs.py", "rdf_partial": "build_rdf_partial_csv_and_figs.py",
         "rdf_sametype": "build_rdf_sametype_csv_and_figs.py", "rdf_cross": "build_rdf_cross_csv_and_figs.py"}

def rep(s, old, new, count=None):
    assert old in s, f"no encontrado: {old[:60]!r}"
    return s.replace(old, new) if count is None else s.replace(old, new, count)

def patch_rdf(s):
    s = rep(s, 'RATES = ["1e12", "1e11", "2e10"]', 'RATES = ["1e12", "1e11", "2e10", "qgrad"]')
    s = rep(s, '"2e10": r"$2\\times10^{10}$ K/s"}', '"2e10": r"$2\\times10^{10}$ K/s", "qgrad": "qgrad (graduada)"}')
    s = rep(s, '"2e10": "#1f3fd6"}', '"2e10": "#1f3fd6", "qgrad": "#7b3fbf"}')
    s = s.replace("plt.subplots(3, 1, figsize=(6.2, 9.0)", "plt.subplots(4, 1, figsize=(6.2, 11.6)")
    s = s.replace('.png"', '_4curvas.png"')
    return s

for kind, script in KINDS.items():
    d = HERE / kind; (d / "raw" / "qgrad").mkdir(parents=True, exist_ok=True)
    for r in ("2e10", "1e11", "1e12"):
        link = d / "raw" / r
        if not link.exists(): link.symlink_to(E4 / kind / "raw" / r)
    q = d / "raw" / "qgrad"
    for nn in range(1, 6):
        for ext in ("rdf.txt", "log"):
            l = q / f"{nn:02d}.{ext}"
            if not l.is_symlink() and not l.exists(): l.symlink_to(E4 / kind / "raw" / "2e10" / f"{nn:02d}.{ext}")
    lab = "17r" if (q / "17r.rdf.txt").exists() else "17u"
    for ext in ("rdf.txt", "log"):
        l = q / f"17.{ext}"
        if l.is_symlink(): l.unlink()
        if (q / f"{lab}.{ext}").exists(): l.symlink_to(f"{lab}.{ext}")
    src = (E4 / kind / script).read_text()
    out = patch_rdf(src)
    if kind == "rdf_cross":   # fig M (emergencia del hombro Zr-Cu, era solo 2e10): se hace para qgrad
        out = rep(out, 'data[("2e10", T)]', 'data[("qgrad", T)]')
        out = rep(out, r'tasa $2\times10^{10}$ K/s', 'perfil qgrad (1500-700 K graduada)')
    (d / script.replace(".py", "_4curvas.py")).write_text(out)
    print(f"{kind}: raw listo (300 K = {lab}), {script.replace('.py', '_4curvas.py')} escrito")

# ---- fig4_style (Fig 4 Mendez): RATE=qgrad
d = HERE / "fig4_style"; d.mkdir(exist_ok=True)
s = (E4 / "fig4_style" / "build_fig4_style.py").read_text()
s = rep(s, 'RATE = "2e10"', 'RATE = "qgrad"')
s = rep(s, r'tasa $2\\times10^{{10}}$ K/s', 'perfil qgrad')
s = s.replace('.png"', '_4curvas.png"')
(d / "build_fig4_style_4curvas.py").write_text(s)

# ---- BADF
d = HERE / "badf"; d.mkdir(exist_ok=True)
s = (E4 / "badf" / "badf.py").read_text()
s = rep(s, 'RATES = ["1e12", "1e11", "2e10"]', 'RATES = ["1e12", "1e11", "2e10", "qgrad"]')
s = rep(s, 'RATE_LABEL = {"1e12": r"$10^{12}$ K/s", "1e11": r"$10^{11}$ K/s", "2e10": r"$2\\times10^{10}$ K/s"}',
        'RATE_LABEL = {"1e12": r"$10^{12}$ K/s", "1e11": r"$10^{11}$ K/s", "2e10": r"$2\\times10^{10}$ K/s", "qgrad": "qgrad"}')
s = rep(s, '"2e10": "#1f3fd6"}', '"2e10": "#1f3fd6", "qgrad": "#7b3fbf"}')
s = rep(s, 'types, pos, L = read_dump(SA / rate / "snaps" / f"{nn:02d}.dump")', 'types, pos, L = read_dump(snap_path(rate, nn))')
s = rep(s, "# ------------------------------------------------------------------ main",
'''def snap_path(rate, nn):
    """originales: structural_analysis/<tasa>/snaps/NN.dump ; qgrad: qgrad/snaps/<label>.dump (17 -> 17r si existe, si no 17u)."""
    if rate != "qgrad":
        return SA.parent / rate / "snaps" / f"{nn:02d}.dump"
    lab = f"{nn:02d}" if nn != 17 else ("17r" if (SA / "snaps" / "17r.dump").exists() else "17u")
    return SA / "snaps" / f"{lab}.dump"


# ------------------------------------------------------------------ main''')
s = rep(s, 'liquid 1200 K (3 rates)', 'liquid 1200 K (4 curves)')
s = rep(s, '''refs = {"ico"''', '''# --- cutoff UNICO compartido (pedido explicito para el numero comparado): media de los rc de vidrio de las 4 curvas.
# Cada curva ya tiene su propio primer minimo (arriba); como el minimo del vidrio es plano (+-0.1 A) y <CN> cambia ~0.2, el
# <CN> y P(60 deg) que se COMPARAN entre curvas se recalculan con este rc unico.
RC_SHARED = {st: float(np.mean([float(m[6]) for m in meta_rows if m[1] == st and m[2] == "all"])) for st in STATES}
shared_rows = []
for rate in RATES:
    for state, (nn, T) in STATES.items():
        types, pos, L = read_dump(snap_path(rate, nn))
        P, meta = badf_config(types, pos, L, RC_SHARED[state])
        for c in ("all", "Cu"):
            m60 = (ALPHA > 55) & (ALPHA < 66)
            shared_rows.append([rate, state, c, f"{RC_SHARED[state]:.3f}", f"{meta[c]['mean_CN']:.3f}",
                                f"{P[c][m60].max():.4f}", f"{ALPHA[m60][np.argmax(P[c][m60])]:.1f}"])
            print(f"[rc unico {RC_SHARED[state]:.3f}] {rate:5s} {state:9s} {c:3s} <CN>={meta[c]['mean_CN']:.3f} P(60)max={P[c][m60].max():.3f}")
with open(HERE / "badf_sharedcutoff.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "state", "center", "rc_shared", "mean_CN", "P_max_55_66deg", "alpha_at_max"])
    w.writerows(shared_rows)

refs = {"ico"''')
s = rep(s, 'ETAPA4 / "fig_C_badf.png"', 'ETAPA4 / "fig_C_badf_4curvas.png"')
(d / "badf_4curvas.py").write_text(s)
(d / "run_badf_4curvas.sbatch").write_text("""#!/bin/bash
# BADF 4 curvas (Python, 1 core). Uso: sbatch run_badf_4curvas.sbatch
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -J qg_badf
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PYTHONDONTWRITEBYTECODE=1
/scratch/negishi/aaabt/venv_structural/bin/python badf_4curvas.py
""")
print("fig4_style y badf listos")
