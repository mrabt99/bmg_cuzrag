#!/usr/bin/env python3
"""Etapa 1b: Tg por interseccion de dos rectas (rama liquida / rama
vitrea) en V(T) y en Ep(T), adaptado de la logica de vol_ene.py del
director pero generalizado:

  - Opera sobre la tabla concatenada de la Tarea B (thermo_<tasa>.csv:
    step,T,Ep,Lx,rho,P), no sobre un log.lammps de una sola corrida.
  - nparticles ya no hace falta (Ep del CSV ya viene por atomo, ver
    Tarea B).
  - En vez del running_mean(50) + downsampling por modulo-60 de
    vol_ene.py (afinado a mano para esa corrida de 8192 atomos),
    binea T en intervalos de ancho fijo (T_BIN) y promedia V y Ep
    dentro de cada bin. Funciona igual sin importar la densidad de
    muestreo original (2e10 muestrea cada 5000 pasos, 1e12 cada 100).
  - Las ventanas de ajuste liquido/vitreo NO estan hardcodeadas por
    indice de array ni copiadas de la corrida vieja: se pasan como
    rangos de T (T_LIQ_MIN, T_LIQ_MAX, T_GLASS_MIN, T_GLASS_MAX),
    elegidos mirando la curva de esta corrida (ver plots
    diagnostico_*.png).
  - Tg se calcula explicitamente como interseccion algebraica de las
    dos rectas ajustadas, no aproximado a ojo con un zoom.
"""
import csv
import numpy as np


def load_thermo(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    T = np.array([float(r["T"]) for r in rows])
    Ep = np.array([float(r["Ep"]) for r in rows])
    Lx = np.array([float(r["Lx"]) for r in rows])
    V = Lx ** 3
    order = np.argsort(-T)  # de caliente a frio
    return T[order], Ep[order], V[order]


def bin_by_temperature(T, y, t_bin):
    """Promedia y en bines de ancho t_bin sobre el rango de T. Devuelve
    (T_binned, y_binned) ordenados de mayor a menor T."""
    t_lo, t_hi = T.min(), T.max()
    edges = np.arange(np.floor(t_lo / t_bin) * t_bin,
                       np.ceil(t_hi / t_bin) * t_bin + t_bin, t_bin)
    T_out, y_out = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (T >= lo) & (T < hi)
        if mask.sum() == 0:
            continue
        T_out.append(T[mask].mean())
        y_out.append(y[mask].mean())
    T_out = np.array(T_out)
    y_out = np.array(y_out)
    order = np.argsort(-T_out)
    return T_out[order], y_out[order]


def fit_line(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    return m, c


TG_WINDOWS = {
    # Elegidas mirando diagnostico_zoom_*.png (kink visible ~600-900K en
    # ambas tasas) y verificando residuos de cada ventana (ver
    # Tg_debug_*.png): la rama liquida NO se extiende hasta 2000K porque
    # Ep(T) tiene una curvatura real (no ruido) en todo ese rango que
    # sesga la pendiente y hunde Tg(Ep) ~200K por debajo de Tg(V) de
    # forma espuria. 900-1200K es la ventana mas angosta pegada a la
    # transicion donde los residuos de Ep(T) ya son ruido chico, no
    # curvatura sistematica.
    "liquid": (900.0, 1200.0),
    "glass": (300.0, 600.0),
}


def two_branch_tg(T_binned, y_binned, t_liq_range, t_glass_range):
    """Ajusta una recta en la rama liquida (t_liq_range) y otra en la
    vitrea (t_glass_range), devuelve (Tg, m_liq, c_liq, m_glass, c_glass)."""
    liq_mask = (T_binned >= t_liq_range[0]) & (T_binned <= t_liq_range[1])
    glass_mask = (T_binned >= t_glass_range[0]) & (T_binned <= t_glass_range[1])
    if liq_mask.sum() < 2 or glass_mask.sum() < 2:
        raise ValueError(
            f"muy pocos puntos para ajustar: liquido={liq_mask.sum()} "
            f"vitreo={glass_mask.sum()} (revisar ventanas T_LIQ/T_GLASS)")
    m_liq, c_liq = fit_line(T_binned[liq_mask], y_binned[liq_mask])
    m_glass, c_glass = fit_line(T_binned[glass_mask], y_binned[glass_mask])
    Tg = (c_glass - c_liq) / (m_liq - m_glass)
    return Tg, (m_liq, c_liq), (m_glass, c_glass)
