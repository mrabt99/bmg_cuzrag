"""Estilo único para TODAS las figuras del proyecto bmg_cuzrag.

Uso mínimo (desde cualquier script, esté donde esté dentro del repo):

    import sys, pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[N] / "estilo"))
    import bmg_estilo as be

    fig, ax = be.figura("media")                   # UN eje, tamaño fijo (aplica rcParams)
    ax.plot(T, y, **be.linea("2e10"))              # color/marcador/etiqueta fijos
    ax.set_xlabel(be.ETQ["T"]); ax.set_ylabel(be.ETQ["Ep"])
    be.guardar(fig, "mq_tg", script=__file__, fuentes=["structural_analysis/2e10/thermo_2e10.csv"])

`guardar` aplica la coma decimal a todos los ejes, revisa la figura (avisa si hay
números con punto decimal en textos, ejes sin etiqueta, títulos de figura, etc.),
escribe SOLO PNG (300 dpi) en /scratch/negishi/aaabt/bitacora/figs/
y registra la figura en bitacora/figs/indice_figuras.csv.

Las reglas y el porqué de cada decisión están en estilo/README.md y en la skill
/scratch/negishi/aaabt/.claude/skills/figuras-bmg/SKILL.md. Si cambiás algo acá, actualizá esos dos.
"""
from __future__ import annotations

import csv
import datetime as _dt
import re
import subprocess
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colors as mcolors  # noqa: E402
from matplotlib import ticker as mticker  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]          # .../aaabt/bmg_cuzrag (repo)
BASE = RAIZ.parent                                  # .../aaabt
# La bitácora vive FUERA del repo (usuario, 2026-11-13: que bmg_cuzrag no pese al bajarlo).
DIR_FIGS = BASE / "bitacora" / "figs"
INDICE = DIR_FIGS / "indice_figuras.csv"

CM = 1 / 2.54

# ---------------------------------------------------------------------------
# Tamaños. Beamer 16:9 = 160 x 90 mm; con título de frame queda ~150 x 68 mm.
# Las figuras se incluyen SIEMPRE a su tamaño natural (width=\figancho{...} en
# el .tex), así la letra de la figura mide lo mismo en todas las slides.
# ---------------------------------------------------------------------------
#
# REGLA: una figura = un eje = un archivo. Nada de subplots: si hacen falta dos
# paneles, son dos archivos y se ponen lado a lado en LaTeX (subcaption pone (a), (b)).
# Por eso casi todo debería ser "media" (dos por fila en la slide y en el informe).
TAMANOS = {
    "media": (7.2, 5.4),         # default: dos lado a lado, o figura + texto
    "dos_tercios": (9.6, 6.0),   # figura + columna angosta de texto
    "ancha": (14.5, 6.0),        # una sola figura a todo el ancho (series largas)
    "cuadrada": (6.0, 6.0),      # mapas, heatmaps
}

# ---------------------------------------------------------------------------
# Paleta. Okabe-Ito (segura para daltonismo). UNA muestra = UN color, en todas
# las figuras del proyecto. Nunca usar el ciclo por defecto de matplotlib.
# ---------------------------------------------------------------------------
MUESTRAS = {
    "2e10": dict(color="#0072B2", marker="o", ls="-",
                 etiqueta=r"$2\times10^{10}$ K/s"),
    "1e11": dict(color="#009E73", marker="s", ls="-",
                 etiqueta=r"$10^{11}$ K/s"),
    "1e12": dict(color="#E69F00", marker="^", ls="-",
                 etiqueta=r"$10^{12}$ K/s"),
    # Perfiles graduados: la etiqueta es corta a propósito (una leyenda más ancha que
    # el panel hace colapsar constrained_layout). El perfil de tasas se explica en la slide.
    "qgrad": dict(color="#CC79A7", marker="D", ls="--", etiqueta="qgrad"),
    "qgrad_ext": dict(color="#D55E00", marker="v", ls="--", etiqueta="qgrad_ext"),
}

# Tasas del barrido 2000->1900 K que no son muestras de producción.
OTRAS_TASAS = {
    "1e10": dict(color="#56B4E9", marker="o", ls=":", etiqueta=r"$10^{10}$ K/s"),
    "5e10": dict(color="#7F7F7F", marker="o", ls=":", etiqueta=r"$5\times10^{10}$ K/s"),
    "2e11": dict(color="#9A9A9A", marker="s", ls=":", etiqueta=r"$2\times10^{11}$ K/s"),
    "5e11": dict(color="#B5B5B5", marker="s", ls=":", etiqueta=r"$5\times10^{11}$ K/s"),
}

# Especies. Se usan SOLO cuando la especie es lo único que cambia en el eje;
# si en el mismo eje varían muestra y especie: color = muestra, marcador/panel = especie.
ESPECIES = {
    "Zr": dict(color="#4D4D4D", marker="o", etiqueta="Zr"),
    "Cu": dict(color="#B5651D", marker="s", etiqueta="Cu"),
    "Ag": dict(color="#6F8FA6", marker="^", etiqueta="Ag"),
}

# Referencias externas: siempre en gris/negro y con marcador hueco.
REFERENCIAS = {
    "mendez": dict(color="0.15", marker="s", mfc="none", ls="none",
                   etiqueta="Méndez y Jaroszewicz"),
    "ryltsev": dict(color="0.40", marker="o", mfc="none", ls="none",
                    etiqueta="Ryltsev et al."),
    "experimento": dict(color="k", marker="*", mfc="k", ls="none",
                        etiqueta="experimental"),
    # Conteos bibliográficos (OpenAlex), figura intro_publicaciones (17/11/2026).
    "biblio_mg": dict(color="k", marker="none", ls="-",
                      etiqueta="«metallic glass»"),
    "biblio_bmg": dict(color="0.55", marker="none", ls="-",
                       etiqueta="«bulk metallic glass»"),
}

# Banda de Tg de referencia (Méndez y Jaroszewicz: 664 K Wendt-Abraham, 692 K Ep).
TG_REF_BANDA = (664.0, 692.0)
# Tg propias (structural_analysis/Tg_resultados.csv, promedio V y Ep). Para marcar en figuras.
# 17/11/2026: 1e12 = corrida unificada (rate_1e12/, structural_analysis/1e12u/Tg_1e12u.csv); la vieja daba 700.6.
TG_PROPIAS = {"2e10": 641.7, "1e11": 666.2, "1e12": 698.0}

# Colormap para curvas coloreadas por temperatura: frío = oscuro, caliente = claro.
# inferno recortado (sin el amarillo pálido final, que no se ve sobre blanco).
CMAP_T = mcolors.LinearSegmentedColormap.from_list(
    "T_bmg", plt.get_cmap("inferno")(__import__("numpy").linspace(0.05, 0.85, 256)))
# Colormap para mapas de densidad / heatmaps de una magnitud >= 0.
CMAP_SEC = "viridis"
# Colormap divergente (diferencias con signo): centrado en 0 con TwoSlopeNorm.
CMAP_DIV = "RdBu_r"

# ---------------------------------------------------------------------------
# Etiquetas de ejes: símbolo en itálica (mathtext), unidades entre paréntesis,
# en castellano. Usar estas en vez de escribirlas a mano.
# ---------------------------------------------------------------------------
ETQ = {
    "T": r"$T$ (K)",
    "t_ns": r"$t$ (ns)",
    "t_ps": r"$t$ (ps)",
    "Ep": r"$E_p$ (eV/átomo)",
    "Ep_meV": r"$E_p$ (meV/átomo)",
    "V_at": r"$V/N$ (Å$^3$/átomo)",
    "rho": r"$\rho$ (g/cm$^3$)",
    "P_GPa": r"$P$ (GPa)",
    "P_bar": r"$P$ (bar)",
    "r": r"$r$ (Å)",
    "g": r"$g(r)$",
    "eps": r"$\varepsilon$",
    "sigma": r"$\sigma$ (GPa)",
    "E": r"$E$ (GPa)",
    "G": r"$G$ (GPa)",
    "kappa": r"$\kappa$ (W m$^{-1}$ K$^{-1}$)",
    "eta": r"$\eta$ (mPa s)",
    "D": r"$D$ (m$^2$/s)",
    "tasa": r"$\dot{T}$ (K/s)",
    "frac_pct": "fracción (%)",
    "phi": r"$\varphi(r)$ (eV)",
}


# ---------------------------------------------------------------------------
# Números con coma decimal
# ---------------------------------------------------------------------------
def num(x, nd=1, mathtext=False):
    """Número con coma decimal y sin separador de miles: num(641.72, 1) -> '641,7'.
    Con mathtext=True devuelve '641{,}7' (para usar adentro de $...$)."""
    s = f"{x:.{nd}f}"
    return s.replace(".", "{,}" if mathtext else ",")


def pm(x, err, nd=2):
    """'0,79 ± 0,02'"""
    return f"{num(x, nd)} ± {num(err, nd)}"


class FormateadorComa(mticker.ScalarFormatter):
    """ScalarFormatter normal (mathtext, offset, notación científica) pero con coma
    decimal. No usa locale: es_AR agrupa miles con punto (1000 K -> '1.000')."""

    def __call__(self, x, pos=None):
        return super().__call__(x, pos).replace(".", "{,}")

    def get_offset(self):
        return super().get_offset().replace(".", "{,}")


def _poner_coma(ax):
    for eje in (ax.xaxis, ax.yaxis):
        f = eje.get_major_formatter()
        if type(f) is mticker.ScalarFormatter:
            nf = FormateadorComa(useOffset=f.get_useOffset(), useMathText=True)
            nf.set_powerlimits((-3, 4))
            eje.set_major_formatter(nf)


# ---------------------------------------------------------------------------
# rcParams
# ---------------------------------------------------------------------------
RC = {
    # Times New Roman (pedido del usuario: es la fuente del informe de avance).
    # En Negishi no está instalada: Liberation Serif tiene sus mismas métricas.
    # Fórmulas en STIX (diseñada para acompañar a Times).
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Liberation Serif", "TeX Gyre Termes",
                   "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.labelsize": 9.5,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "legend.handlelength": 1.8,
    "legend.borderaxespad": 0.4,
    "axes.linewidth": 0.8,
    "axes.formatter.use_mathtext": True,
    "axes.formatter.limits": (-3, 4),
    "axes.prop_cycle": matplotlib.cycler(color=["#000000"]),  # obliga a elegir color
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2.0,
    "ytick.minor.size": 2.0,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    "errorbar.capsize": 2,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "figure.dpi": 150,
    "pdf.fonttype": 42,
    "image.cmap": "viridis",
}


def aplicar():
    plt.rcParams.update(RC)


def figura(tamano="media", **kw):
    """(fig, ax) con UN solo eje y uno de los tamaños fijos (cm) de TAMANOS."""
    aplicar()
    w, h = TAMANOS[tamano]
    kw.setdefault("constrained_layout", True)
    return plt.subplots(1, 1, figsize=(w * CM, h * CM), **kw)


def linea(clave, **extra):
    """kwargs para ax.plot / ax.errorbar de una muestra, tasa, especie o referencia."""
    for tabla in (MUESTRAS, OTRAS_TASAS, ESPECIES, REFERENCIAS):
        if clave in tabla:
            d = dict(tabla[clave])
            break
    else:
        raise KeyError(f"'{clave}' no está en la paleta del proyecto (bmg_estilo.py)")
    d["label"] = d.pop("etiqueta")
    d.update(extra)
    return d


def leyenda(ax, *args, afuera=False, **kw):
    """ax.legend que no deforma el layout. Una leyenda interna más ancha que el panel
    hace que constrained_layout achique los ejes hasta colapsarlos; acá se la excluye
    del cálculo (set_in_layout(False)). Si no entra, que se vea en el PNG y se mueva.
    afuera=True: a la derecha del eje (sí participa del layout)."""
    if afuera:
        kw.setdefault("loc", "upper left")
        kw.setdefault("bbox_to_anchor", (1.02, 1.0))
        return ax.legend(*args, **kw)
    leg = ax.legend(*args, **kw)
    leg.set_in_layout(False)
    return leg


def banda_tg_ref(ax, **kw):
    """Banda gris de Tg de Méndez y Jaroszewicz (664–692 K) en un eje con T en x."""
    kw.setdefault("color", "0.85")
    kw.setdefault("zorder", 0)
    kw.setdefault("lw", 0)
    return ax.axvspan(*TG_REF_BANDA, **kw)


def marcar_tg(ax, muestra, tg=None, **kw):
    """Línea vertical punteada en Tg de una muestra, con su color."""
    tg = TG_PROPIAS[muestra] if tg is None else tg
    kw.setdefault("ls", ":")
    kw.setdefault("lw", 1.0)
    return ax.axvline(tg, color=MUESTRAS[muestra]["color"], **kw)


# ---------------------------------------------------------------------------
# Revisión y guardado
# ---------------------------------------------------------------------------
_RE_PUNTO = re.compile(r"(?<![\w{])\d+\.\d+")


def revisar(fig):
    """Devuelve una lista de avisos sobre cosas que rompen el estilo del proyecto."""
    avisos = []
    ejes = [ax for ax in fig.axes if ax.get_label() != "<colorbar>"]
    if len(ejes) > 1:
        avisos.append(f"tiene {len(ejes)} ejes: una figura = un eje = un archivo "
                      "(los paneles se arman en LaTeX)")
    if fig._suptitle is not None and fig._suptitle.get_text():
        avisos.append("tiene suptitle: el título va en el frame de Beamer, no en la figura")
    textos = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible()]
    for t in textos:
        s = t.get_text()
        if _RE_PUNTO.search(s.replace(r"\mathdefault", "")):
            avisos.append(f"número con punto decimal en texto: {s!r} (usar be.num())")
    for i, ax in enumerate(fig.axes):
        if not ax.get_visible() or ax.get_label() == "<colorbar>":
            continue
        if not ax.axison:
            continue
        for nombre, lab in (("x", ax.get_xlabel()), ("y", ax.get_ylabel())):
            if not lab and getattr(ax, f"get_{nombre}ticklabels")() and \
                    any(tl.get_text() for tl in getattr(ax, f"get_{nombre}ticklabels")()):
                avisos.append(f"eje {i}: falta etiqueta en {nombre}")
            elif lab and "(" not in lab and not re.search(r"\$[^$]*\$\s*$", lab) \
                    and "%" not in lab:
                avisos.append(f"eje {i}: etiqueta {nombre} sin unidades entre paréntesis: {lab!r}"
                              " (si es adimensional, ignorar)")
        if ax.get_title():
            avisos.append(f"eje {i}: tiene título ({ax.get_title()!r}); el título va en el"
                          " caption/slide, no en la figura")
    return avisos


def _git_rev():
    try:
        return subprocess.run(["git", "-C", str(RAIZ), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:
        return ""


def guardar(fig, nombre, script, fuentes=(), destino=None, estricto=False):
    """Coma decimal + revisión + PNG (300 dpi, solo PNG: pedido del usuario) + índice.

    nombre : sin extensión, prefijo de bloque (p. ej. 'mq_tg', 'est_voronoi_cu').
    script : pasar __file__.
    fuentes: rutas (relativas a bmg_cuzrag/) de los datos que lee la figura.
    destino: carpeta de salida; por defecto bitacora/figs/ (y se registra en su índice).
             Un script de análisis dentro del repo puede pasar su propia carpeta.
    """
    for ax in fig.axes:
        _poner_coma(ax)
    avisos = revisar(fig)
    for a in avisos:
        warnings.warn(f"[{nombre}] {a}", stacklevel=2)
    if estricto and avisos:
        raise RuntimeError(f"{nombre}: {len(avisos)} avisos de estilo")

    dest = DIR_FIGS if destino is None else Path(destino)
    dest.mkdir(parents=True, exist_ok=True)
    try:
        script_rel = str(Path(script).resolve().relative_to(BASE))
    except ValueError:
        script_rel = str(Path(script).resolve())
    meta = {"Software": f"{script_rel} | datos: {'; '.join(fuentes)}"}
    fig.savefig(dest / f"{nombre}.png", metadata=meta)
    if dest != DIR_FIGS:
        print(f"[bmg_estilo] {dest / nombre}.png  ({len(avisos)} avisos)")
        return avisos

    filas = {}
    if INDICE.exists():
        with open(INDICE, newline="") as f:
            filas = {r["figura"]: r for r in csv.DictReader(f)}
    filas[nombre] = {
        "figura": nombre,
        "script": script_rel,
        "fuentes": "; ".join(fuentes),
        "fecha": _dt.date.today().isoformat(),
        "git": _git_rev(),
        "avisos": len(avisos),
    }
    with open(INDICE, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(next(iter(filas.values())).keys()))
        w.writeheader()
        for k in sorted(filas):
            w.writerow(filas[k])
    print(f"[bmg_estilo] {nombre}.png  ({len(avisos)} avisos)")
    return avisos
