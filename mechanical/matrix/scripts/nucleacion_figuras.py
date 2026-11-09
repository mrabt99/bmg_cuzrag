#!/usr/bin/env python3
"""Figuras/resumen del mapeo de nucleacion (lee nucleacion_casos.csv y nucleacion_indices.csv; escribe figs/nucleacion_*.png y nucleacion_resumen.txt). Uso (venv): python nucleacion_figuras.py"""
import os, csv
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); F = f"{M}/figs"
C = list(csv.DictReader(open(f"{M}/nucleacion_casos.csv"))); I = list(csv.DictReader(open(f"{M}/nucleacion_indices.csv")))
COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad_ext": "#b5338a"}; S = ["2e10", "1e11", "1e12", "qgrad_ext"]
MOD = ["AUC_M0_especie", "AUC_M1_+CN", "AUC_M2_+ico", "AUC_M3_+vol", "AUC_M4_+entorno"]; ML = ["especie", "+CN", "+ico", "+vol. Voronoi", "+entorno 6 Å"]
def clean(a):
    a.grid(color="#e6e5e1", lw=.7); a.set_axisbelow(True)
    for s in ("top", "right"): a.spines[s].set_visible(False)
def agg(s, T, er, k): v = [float(r[k]) for r in C if r["sample"] == s and float(r["T"]) == T and r["erate"] == er]; return (np.mean(v), np.std(v, ddof=1)) if len(v) > 1 else (np.nan, np.nan)
L = []
fg, ax = plt.subplots(2, 3, figsize=(17, 9.5), constrained_layout=True)
# (0,0) AUC por modelo a 300K 1e9
for i, s in enumerate(S):
    ax[0, 0].errorbar(np.arange(5) + (i - 1.5) * .12, [agg(s, 300, "1e9", m)[0] for m in MOD], yerr=[agg(s, 300, "1e9", m)[1] for m in MOD], color=COL[s], marker="o", capsize=2, label=s)
ax[0, 0].axhline(.5, color="#8a8985", lw=.8); ax[0, 0].set_xticks(range(5)); ax[0, 0].set_xticklabels(ML, fontsize=8.5); ax[0, 0].set(ylabel="AUC (entrena x < L/2, prueba x ≥ L/2)", title="¿La estructura inicial predice dónde nace la actividad? (300 K, 1e9)"); ax[0, 0].legend(frameon=False, fontsize=8.5)
# (0,1) coeficientes M4
cn = ["Cu", "Ag", "CN", "ico", "vol_rel", "cg_ico", "cg_vol"]
for i, s in enumerate(S): ax[0, 1].bar(np.arange(7) + (i - 1.5) * .2, [agg(s, 300, "1e9", "coef_" + k)[0] for k in cn], .2, yerr=[agg(s, 300, "1e9", "coef_" + k)[1] for k in cn], color=COL[s], error_kw=dict(lw=.7, capsize=1.5))
ax[0, 1].axhline(0, color="#8a8985", lw=.8); ax[0, 1].set_xticks(range(7)); ax[0, 1].set_xticklabels(["Cu", "Ag", "CN", "ico", "vol.", "ico en\nr<6 Å", "vol. en\nr<6 Å"], fontsize=8.5); ax[0, 1].set(ylabel="coeficiente logístico (variables estandarizadas)", title="Modelo completo: qué empuja a un átomo a nuclear (300 K, 1e9)")
# (0,2) correlacion voxel
for i, s in enumerate(S):
    for j, (er, mk) in enumerate((("1e8", "s"), ("1e9", "o"), ("1e10", "^"))):
        m, sd = agg(s, 300, er, "rho_vox_ico_act"); ax[0, 2].errorbar(i + (j - 1) * .15, m, yerr=sd, color=COL[s], marker=mk, capsize=2, mfc=COL[s] if er == "1e9" else "none")
ax[0, 2].axhline(0, color="#8a8985", lw=.8); ax[0, 2].set_xticks(range(4)); ax[0, 2].set_xticklabels(S); ax[0, 2].set(ylabel="Spearman (fracción ico inicial vs fracción activa), voxeles ~10 Å", title="Correlación espacial (cuadrado/círculo/triángulo = 1e8/1e9/1e10)")
# (1,0) AUC vs T
for s in ["2e10", "1e11", "1e12"]:
    Ts = sorted({float(r["T"]) for r in C if r["sample"] == s and r["erate"] == "1e9"}); ax[1, 0].errorbar(Ts, [agg(s, T, "1e9", MOD[4])[0] for T in Ts], yerr=[agg(s, T, "1e9", MOD[4])[1] for T in Ts], color=COL[s], marker="o", capsize=2, label=s)
ax[1, 0].set(xlabel="T (K)", ylabel="AUC del modelo completo", title="Poder predictivo vs temperatura (1e9)"); ax[1, 0].axhline(.5, color="#8a8985", lw=.8); ax[1, 0].legend(frameon=False, fontsize=8.5)
# (1,1) coef ico vs T
for s in ["2e10", "1e11", "1e12"]:
    Ts = sorted({float(r["T"]) for r in C if r["sample"] == s and r["erate"] == "1e9"}); ax[1, 1].errorbar(Ts, [agg(s, T, "1e9", "coef_ico")[0] for T in Ts], yerr=[agg(s, T, "1e9", "coef_ico")[1] for T in Ts], color=COL[s], marker="o", capsize=2, label=s)
ax[1, 1].axhline(0, color="#8a8985", lw=.8); ax[1, 1].set(xlabel="T (K)", ylabel="coeficiente de 'es ico' (controlando CN, volumen, entorno)", title="Efecto protector del icosaedro vs T (1e9)")
# (1,2) riesgo relativo por indice, Cu y Zr, 300K 1e9 agrupado
agg_i = defaultdict(lambda: defaultdict(lambda: [0, 0.0]))
for r in I:
    if float(r["T"]) == 300 and r["erate"] == "1e9":
        a = agg_i[(r["especie"], r["indice"])][r["sample"]]; a[0] += int(r["N"]); a[1] += float(r["riesgo_relativo"]) * int(r["N"])
top = sorted(agg_i, key=lambda k: -sum(v[0] for v in agg_i[k].values()))
sel = [k for k in top if k[0] == "Cu"][:7] + [k for k in top if k[0] == "Zr"][:7]
for i, s in enumerate(S): ax[1, 2].barh(np.arange(len(sel)) + (i - 1.5) * .2, [agg_i[k][s][1] / agg_i[k][s][0] if agg_i[k][s][0] else np.nan for k in sel], .2, color=COL[s])
ax[1, 2].axvline(1, color="#8a8985", lw=.8); ax[1, 2].set_yticks(range(len(sel))); ax[1, 2].set_yticklabels([f"{k[0]} <{k[1]}>" for k in sel], fontsize=8.5); ax[1, 2].invert_yaxis(); ax[1, 2].set(xlabel="riesgo relativo de nuclear (vs su especie)", title="Poliedros de Voronoi más frecuentes (300 K, 1e9)")
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
hs = [Patch(color=COL[s], label=s) for s in S]
ax[0, 1].legend(handles=hs, frameon=False, fontsize=8.5, title="vidrio"); ax[1, 2].legend(handles=hs, frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(1.01, 1.0), title="vidrio")
ax[0, 2].legend(handles=[Line2D([], [], color="k", marker=m, ls="", mfc="none", label=l) for m, l in (("s", "1e8"), ("o", "1e9"), ("^", "1e10"))], title="tasa de deformación (s⁻¹)", frameon=False, fontsize=8.5, loc="lower left")
ax[1, 1].legend(frameon=False, fontsize=8.5, title="vidrio"); ax[0, 0].legend(frameon=False, fontsize=8.5, title="vidrio")
for a in ax.flat: a.title.set_fontsize(9.5)
for a in ax.flat: clean(a)
fg.suptitle("Mapeo local de nucleación: estructura Voronoi en ε = 0 vs actividad no afín (D²min) hasta el pico", fontsize=11.5); fg.savefig(f"{F}/nucleacion_resumen.png", dpi=135); plt.close(fg)
# numeros
L.append("== AUC (media x,y,z) 300K 1e9: " + " | ".join(f"{s}: " + " ".join(f"{agg(s,300,'1e9',m)[0]:.3f}" for m in MOD) for s in S))
L.append("== fraccion activa antes del pico: " + " ".join(f"{s} {agg(s,300,'1e9','frac_activo')[0]*100:.1f}%" for s in S))
L.append("== coef M4 (300K 1e9): " + " | ".join(f"{s}: " + " ".join(f"{k}={agg(s,300,'1e9','coef_'+k)[0]:+.2f}" for k in cn) for s in S))
L.append("== Spearman voxel ico-actividad (1e9): " + " ".join(f"{s} {agg(s,300,'1e9','rho_vox_ico_act')[0]:+.2f}" for s in S))
L.append("== AUC M4 vs T (1e9): " + " | ".join(f"{s}: " + " ".join(f"{int(T)}K={agg(s,T,'1e9',MOD[4])[0]:.3f}" for T in sorted({float(r['T']) for r in C if r['sample']==s and r['erate']=='1e9'})) for s in ["2e10","1e11","1e12"]))
L.append("== coef ico vs T (1e9): " + " | ".join(f"{s}: " + " ".join(f"{int(T)}K={agg(s,T,'1e9','coef_ico')[0]:+.2f}" for T in sorted({float(r['T']) for r in C if r['sample']==s and r['erate']=='1e9'})) for s in ["2e10","1e11","1e12"]))
L.append("== riesgo relativo por indice (300K,1e9; ponderado por N): " + "; ".join(f"{k[0]}<{k[1]}> " + "/".join(f"{agg_i[k][s][1]/agg_i[k][s][0]:.2f}" if agg_i[k][s][0] else "-" for s in S) for k in sel))
open(f"{M}/nucleacion_resumen.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
