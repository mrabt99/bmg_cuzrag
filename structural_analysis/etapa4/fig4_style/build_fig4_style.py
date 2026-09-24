#!/usr/bin/env python3
"""Etapa 4 / Tarea 6 -- replica del estilo de Figura4_Mendez.png (../../../referencias/) para
los pares Cu-Cu, Zr-Zr y Zr-Cu (cruzado). Panel izquierdo: g(r)+offset apiladas por T (waterfall).
Panel derecho: zoom del 2do pico, sin offset, para ver la evolucion/split con T.

Entrada : ../rdf_sametype/rdf_sametype_all.csv (g_CuCu, g_ZrZr) y ../rdf_cross/rdf_cross_all.csv (g_ZrCu)
Salida  : ../fig_N_fig4style_CuCu.png, ../fig_O_fig4style_ZrZr.png, ../fig_P_fig4style_ZrCu.png
Tasa: 2e10 (misma que fig_J / fig_M). Sin pandas: numpy + csv.
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
RATE = "2e10"


def load(path, col):
    data = {}
    with open(path) as f:
        rd = csv.DictReader(f)
        for row in rd:
            if row["rate"] != RATE:
                continue
            T = int(float(row["T"]))
            data.setdefault(T, {"r": [], "g": []})
            data[T]["r"].append(float(row["r"]))
            data[T]["g"].append(float(row[col]))
    for T in data:
        data[T]["r"] = np.array(data[T]["r"])
        data[T]["g"] = np.array(data[T]["g"])
    return data


DATASETS = [
    ("CuCu", HERE.parent / "rdf_sametype" / "rdf_sametype_all.csv", "g_CuCu",
     r"$g_{\mathrm{Cu-Cu}}(r)$", "fig_N_fig4style_CuCu.png"),
    ("ZrZr", HERE.parent / "rdf_sametype" / "rdf_sametype_all.csv", "g_ZrZr",
     r"$g_{\mathrm{Zr-Zr}}(r)$", "fig_O_fig4style_ZrZr.png"),
    ("ZrCu", HERE.parent / "rdf_cross" / "rdf_cross_all.csv", "g_ZrCu",
     r"$g_{\mathrm{Zr-Cu}}(r)$", "fig_P_fig4style_ZrCu.png"),
]

OFFSET_STEP = 0.45  # separacion vertical entre curvas sucesivas (panel izq.)
XMAX = 9.9


CMAP = plt.get_cmap("nipy_spectral")


def make_colors(Ts):
    """Rainbow continuo (alta T = violeta/azul, baja T = rojo), igual sentido que Fig 4 Mendez
    (300 K queda arriba y con tonos calidos, 2000/1900 K abajo con tonos frios)."""
    Tmin, Tmax = min(Ts), max(Ts)
    # invertido: T baja -> extremo "caliente" del cmap (mayor valor de x), T alta -> extremo frio
    return {T: CMAP(0.08 + 0.88 * (1 - (T - Tmin) / (Tmax - Tmin))) for T in Ts}


for label, path, col, ylabel, fname in DATASETS:
    data = load(path, col)
    Ts = sorted(data.keys())  # 300 .. 1900, ascendente
    colors = make_colors(Ts)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 6.5))

    # ---- panel izquierdo: waterfall con offset, T alta abajo (offset 0) -> T baja arriba ----
    Ts_desc_offset = sorted(Ts, reverse=True)  # 1900 primero -> offset 0
    for i, T in enumerate(Ts_desc_offset):
        r, g = data[T]["r"], data[T]["g"]
        m = r <= XMAX
        axL.plot(r[m], g[m] + i * OFFSET_STEP, color=colors[T], lw=1.0)
        y_end = g[m][-1] + i * OFFSET_STEP
        axL.text(XMAX + 0.15, y_end, f"{T} K", color=colors[T], fontsize=7, va="center")
    axL.set_xlim(0, XMAX + 1.3)
    axL.set_xlabel(r"$r$, $\mathrm{\AA}$")
    axL.set_ylabel(ylabel + r" $+$ offset")
    axL.set_title(f"{ylabel} a distintas T (offset {OFFSET_STEP}/curva)", fontsize=10)

    # ---- panel derecho: zoom 2do pico, sin offset, con colorbar (las curvas convergen a la
    # derecha, etiquetar cada una ahi se superpone) ----
    for T in Ts_desc_offset:
        r, g = data[T]["r"], data[T]["g"]
        m = (r >= 3.5) & (r <= XMAX)
        axR.plot(r[m], g[m], color=colors[T], lw=1.0)
    axR.set_xlim(3.5, XMAX)
    axR.set_xlabel(r"$r$, $\mathrm{\AA}$")
    axR.set_ylabel(ylabel)
    axR.set_title("Zoom 2do pico (sin offset)", fontsize=10)
    # colorbar con la MISMA funcion T->color que make_colors (T baja=calido/rojo, T alta=frio/violeta)
    Tmin, Tmax = min(Ts), max(Ts)
    def colorval(T):
        return 0.08 + 0.88 * (1 - (T - Tmin) / (Tmax - Tmin))
    sm = ScalarMappable(cmap=CMAP, norm=Normalize(vmin=0, vmax=1))
    cb = fig.colorbar(sm, ax=axR)
    tick_Ts = [T for T in [300, 600, 900, 1200, 1500, 1900] if T in data]
    cb.set_ticks([colorval(T) for T in tick_Ts])
    cb.set_ticklabels([str(T) for T in tick_Ts])
    cb.set_label(r"$T$, K")

    fig.suptitle(f"Estilo Figura 4 Mendez -- {label}, tasa $2\\times10^{{10}}$ K/s", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / fname, dpi=200)
    plt.close(fig)
    print(f"{fname} escrita")
