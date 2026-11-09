#!/usr/bin/env python3
"""Etapa 1 para qgrad: Tg del PERFIL graduado (doble ajuste lineal, misma logica y ventanas que
compute_tg.py: liquido 900-1200 K, vitreo 300-600 K, bines de 10 K) + comparativo 4 curvas.
NO usa run_tg_1e11.py (sobreescribe Tg_resultados.csv/Tg_comparativo.png); todo se escribe en qgrad/thermo/.
Ejecutar con PYTHONDONTWRITEBYTECODE=1 desde qgrad/thermo/.

AVISO: qgrad no tiene una tasa unica (2e10->1.5e9 K/s entre 1500 y 700 K, 1e11 K/s entre 700 y 300 K).
Su 'Tg' es la del perfil: el codo de V(T)/Ep(T) cae cerca del salto de tasa a 700 K, asi que mezcla
la historia lenta con la rapida. Se reporta tambien la sensibilidad a las ventanas.
"""
import sys, csv
sys.path.insert(0, "../..")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from compute_tg import load_thermo, bin_by_temperature, two_branch_tg, TG_WINDOWS

T_BIN = 10.0
CSV = {"2e10": "../../2e10/thermo_2e10.csv", "1e11": "../../1e11/thermo_1e11.csv",
       "1e12": "../../1e12/thermo_1e12.csv", "qgrad": "thermo_qgrad.csv"}
COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad": "#7b3fbf"}
ALT = {"ventanas estandar (liq 900-1200 / vit 300-600)": TG_WINDOWS,
       "liq 900-1100 / vit 300-500": {"liquid": (900., 1100.), "glass": (300., 500.)},
       "liq 800-1000 / vit 300-600": {"liquid": (800., 1000.), "glass": (300., 600.)},
       "liq 1000-1200 / vit 300-600": {"liquid": (1000., 1200.), "glass": (300., 600.)}}

data, out = {}, []
for k, p in CSV.items():
    T, Ep, V = load_thermo(p)
    Tb, Vb = bin_by_temperature(T, V, T_BIN); _, Epb = bin_by_temperature(T, Ep, T_BIN)
    data[k] = (Tb, Vb, Epb)
    for name, w in ALT.items():
        try:
            tv = two_branch_tg(Tb, Vb, w["liquid"], w["glass"])[0]; te = two_branch_tg(Tb, Epb, w["liquid"], w["glass"])[0]
        except ValueError as e:
            tv = te = float("nan")
        out.append((k, name, tv, te))
with open("Tg_qgrad_y_sensibilidad.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["curva", "ventanas", "Tg_V_K", "Tg_Ep_K"])
    w.writerows((k, n, f"{a:.1f}", f"{b:.1f}") for k, n, a, b in out)
for k, n, a, b in out:
    print(f"{k:6s} {n:48s} Tg(V)={a:6.1f}  Tg(Ep)={b:6.1f}")

fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
for k in CSV:
    Tb, Vb, Epb = data[k]
    tv, te = [x for kk, n, *x in out if kk == k and n.startswith("ventanas estandar")][0]
    ax[0].plot(Tb, Vb, ".", ms=2.5, color=COL[k], alpha=.7, label=f"{k}: Tg(V)={tv:.0f} K"); ax[0].axvline(tv, color=COL[k], ls="--", lw=1)
    ax[1].plot(Tb, Epb, ".", ms=2.5, color=COL[k], alpha=.7, label=f"{k}: Tg(Ep)={te:.0f} K"); ax[1].axvline(te, color=COL[k], ls="--", lw=1)
Tq, Vq, Eq = data["qgrad"]; T2, V2, E2 = data["2e10"]
common = [t for t in Tq if np.any(np.isclose(T2, t, atol=0.6))]
for k in ["1e11", "1e12", "qgrad"]:
    Tk, _, Ek = data[k]
    d = [(t, e - E2[np.argmin(abs(T2 - t))]) for t, e in zip(Tk, Ek) if t <= 1500 and abs(T2[np.argmin(abs(T2 - t))] - t) < 3]
    ax[2].plot([a for a, _ in d], [1000 * b for _, b in d], "-", color=COL[k], label=f"{k} - 2e10")
ax[2].axhline(0, color="k", lw=.6); ax[2].axvline(700, color="grey", ls=":", lw=1); ax[2].text(705, ax[2].get_ylim()[0], " salto de tasa\n 1.5e9->1e11", fontsize=8, va="bottom")
ax[0].set(title="V(T) (bines 10 K)", xlabel="T (K)", ylabel="V (A^3)"); ax[1].set(title="Ep(T) (bines 10 K)", xlabel="T (K)", ylabel="Ep (eV/at)")
ax[2].set(title="Ep(T) relativo a 2e10, T<=1500 K", xlabel="T (K)", ylabel="Ep - Ep(2e10)  (meV/at)")
for a in ax: a.legend(fontsize=8); a.grid(alpha=.25)
ax[0].set_xlim(250, 1600); ax[1].set_xlim(250, 1600); ax[2].set_xlim(250, 1500)
fig.suptitle("Tg comparativo 4 curvas (qgrad = perfil graduado, sin tasa unica)")
fig.tight_layout(); fig.savefig("Tg_comparativo_4curvas.png", dpi=130)
