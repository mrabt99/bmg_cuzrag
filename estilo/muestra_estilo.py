"""Hojas de muestra del estilo (datos inventados, solo para ver colores y tipografía).
Salidas (una figura por archivo, PNG): bitacora/figs/estilo_muestras.*, estilo_especies.*,
estilo_cmapT.*  — volver a correr después de tocar bmg_estilo.py."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bmg_estilo as be  # noqa: E402

# 1) las 5 muestras + Tg
fig, ax = be.figura("media")
T = np.linspace(300, 1200, 10)
for i, k in enumerate(be.MUESTRAS):
    y = 0.80 + 0.02 * i + 1e-4 * (T - 300) * (1 + 0.2 * i)
    ax.errorbar(T, y, yerr=0.008, **be.linea(k))
be.banda_tg_ref(ax)
be.marcar_tg(ax, "2e10")
ax.set_xlabel(be.ETQ["T"])
ax.set_ylabel(be.ETQ["kappa"])
be.leyenda(ax, loc="upper left")
be.guardar(fig, "estilo_muestras", script=__file__)

# 2) especies + referencias externas
fig, ax = be.figura("media")
x = np.arange(3)
for j, (k, d) in enumerate(be.ESPECIES.items()):
    ax.bar(x + 0.27 * (j - 1), [12 + j, 13 + j, 14 + j], width=0.25, color=d["color"],
           label=d["etiqueta"])
for k in ("mendez", "ryltsev"):
    ax.plot(x, [13, 14.5, 15.5] if k == "mendez" else [12.5, 14, 15], **be.linea(k))
ax.set_xticks(x, ["300 K", "700 K", "1200 K"])
ax.minorticks_off()
ax.set_ylabel("número de coordinación")
ax.set_ylim(10, 20)
be.leyenda(ax, loc="upper left", ncol=2)
be.guardar(fig, "estilo_especies", script=__file__)

# 3) curvas coloreadas por T
fig, ax = be.figura("media")
r = np.linspace(2, 8, 300)
temps = np.arange(300, 2000, 200)
for Tk in temps:
    g = 1 + np.exp(-(r - 3.1) ** 2 / 0.05) * (2.6 - Tk / 1500) * np.exp(-(r - 3.1) / 3)
    ax.plot(r, g, color=be.CMAP_T((Tk - temps.min()) / np.ptp(temps)), lw=1)
sm = be.plt.cm.ScalarMappable(cmap=be.CMAP_T, norm=be.plt.Normalize(temps.min(), temps.max()))
cb = fig.colorbar(sm, ax=ax, pad=0.02)
cb.set_label(be.ETQ["T"])
ax.set_xlabel(be.ETQ["r"])
ax.set_ylabel(be.ETQ["g"])
be.guardar(fig, "estilo_cmapT", script=__file__)
