#!/usr/bin/env python3
"""Etapa 1c: aplica el mismo Tg de dos rectas (compute_tg.py) a la
tasa 1e11, con las mismas ventanas ya elegidas para 2e10/1e12
(liquido 900-1200K, vitreo 300-600K -- ver TG_WINDOWS en
compute_tg.py), y regenera Tg_comparativo.png con las 3 tasas.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compute_tg import load_thermo, bin_by_temperature, two_branch_tg, TG_WINDOWS

T_BIN = 10.0

def diagnostico(rate, csv_path, png_path):
    T, Ep, V = load_thermo(csv_path)
    Tb, Vb = bin_by_temperature(T, V, T_BIN)
    _, Epb = bin_by_temperature(T, Ep, T_BIN)

    Tg_V, *_ = two_branch_tg(Tb, Vb, TG_WINDOWS["liquid"], TG_WINDOWS["glass"])
    Tg_Ep, *_ = two_branch_tg(Tb, Epb, TG_WINDOWS["liquid"], TG_WINDOWS["glass"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(Tb, Vb, ".", ms=3)
    axes[0].set_title(f"{rate}: V(T) binned (10K)")
    axes[0].set_xlabel("T (K)"); axes[0].set_ylabel("V (A^3)")
    axes[1].plot(Tb, Epb, ".", ms=3)
    axes[1].set_title(f"{rate}: Ep(T) binned (10K)")
    axes[1].set_xlabel("T (K)"); axes[1].set_ylabel("Ep (eV/at)")
    fig.tight_layout()
    fig.savefig(png_path, dpi=130)
    plt.close(fig)
    return Tg_V, Tg_Ep, Tb, Vb, Epb

Tg_V_1e11, Tg_Ep_1e11, Tb_1e11, Vb_1e11, Epb_1e11 = diagnostico(
    "1e11", "1e11/thermo_1e11.csv", "1e11/diagnostico_1e11.png")
print(f"1e11: Tg(V)={Tg_V_1e11:.1f}K  Tg(Ep)={Tg_Ep_1e11:.1f}K")

# actualiza Tg_resultados.csv (reemplaza la fila 1e11 si ya existiera)
rows = []
with open("Tg_resultados.csv") as f:
    reader = csv.reader(f)
    header = next(reader)
    for r in reader:
        if r[0] != "1e11":
            rows.append(r)
rows.append(["1e11", f"{Tg_V_1e11:.1f}", f"{Tg_Ep_1e11:.1f}"])
rate_order = {"2e10": 0, "1e11": 1, "1e12": 2}
rows.sort(key=lambda r: rate_order.get(r[0], 99))
with open("Tg_resultados.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(header)
    w.writerows(rows)
print("Tg_resultados.csv actualizado:")
print(open("Tg_resultados.csv").read())

# Tg_comparativo.png con las 3 tasas
RATE_STYLE = {"2e10": ("tab:blue", "-"), "1e11": ("tab:orange", "-"), "1e12": ("tab:red", "-")}
RATE_CSV = {"2e10": "2e10/thermo_2e10.csv", "1e11": "1e11/thermo_1e11.csv", "1e12": "1e12/thermo_1e12.csv"}

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for rate in ["2e10", "1e11", "1e12"]:
    T, Ep, V = load_thermo(RATE_CSV[rate])
    Tb, Vb = bin_by_temperature(T, V, T_BIN)
    _, Epb = bin_by_temperature(T, Ep, T_BIN)
    Tg_V, *_ = two_branch_tg(Tb, Vb, TG_WINDOWS["liquid"], TG_WINDOWS["glass"])
    Tg_Ep, *_ = two_branch_tg(Tb, Epb, TG_WINDOWS["liquid"], TG_WINDOWS["glass"])
    color, ls = RATE_STYLE[rate]
    axes[0].plot(Tb, Vb, ".", ms=2, color=color, alpha=0.6)
    axes[0].axvline(Tg_V, color=color, ls="--", lw=1, label=f"{rate}: Tg(V)={Tg_V:.0f}K")
    axes[0].plot(Tg_V, np.interp(Tg_V, Tb[::-1], Vb[::-1]), "x", color=color, ms=10, mew=2)
    axes[1].plot(Tb, Epb, ".", ms=2, color=color, alpha=0.6)
    axes[1].axvline(Tg_Ep, color=color, ls="--", lw=1, label=f"{rate}: Tg(Ep)={Tg_Ep:.0f}K")
    axes[1].plot(Tg_Ep, np.interp(Tg_Ep, Tb[::-1], Epb[::-1]), "x", color=color, ms=10, mew=2)

axes[0].set_title("Tg por volumen"); axes[0].set_xlabel("T (K)"); axes[0].set_ylabel("V (A^3)"); axes[0].legend()
axes[1].set_title("Tg por energia potencial"); axes[1].set_xlabel("T (K)"); axes[1].set_ylabel("Ep (eV/at)"); axes[1].legend()
fig.suptitle("Tg comparativo -- ventanas: liquido 900-1200K, vitreo 300-600K")
fig.tight_layout()
fig.savefig("Tg_comparativo.png", dpi=130)
plt.close(fig)
print("Tg_comparativo.png regenerado con 2e10/1e11/1e12")
