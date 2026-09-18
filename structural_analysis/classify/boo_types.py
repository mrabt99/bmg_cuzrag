#!/usr/bin/env python3
"""Etapa 2b - definiciones compartidas para classify_clusters.

Lee los dumps de BOO ya calculados en ../boo/<tasa>/NN.boo.dump (NO los
regenera) y define los puntos ideales de la Tabla I de Ryltsev et al.
2016 mas los criterios de clasificacion.

Layout de columnas del dump (confirmado en la Tarea E, ver
../boo/validate_boo.lmp): id type + 4 tandas de [Q4, Q6, W4hat, W6hat],
una tanda por cada nnn = 12, 14, 15, 16, en ese orden.

Tipos de atomo LAMMPS: los inputs de produccion usan
`pair_coeff * * ZrCuAg.lammps.eam Zr Cu Ag` -> 1=Zr, 2=Cu, 3=Ag
(mismo mapeo que asume boo.py de Jaroszewicz).

El w4/w6 de la Tabla I es el invariante NORMALIZADO, que corresponde
exactamente a nuestro What4/What6 (columnas de `wl/hat yes`), no al
Wl crudo -- confirmado por el usuario antes de escribir esto.

--- REVISION (Tarea F) -------------------------------------------------
El primer intento clasificaba hcp/fcc/Z14/Z15 por bola de distancia en
(q4, q6, What6) con un epsilon unico. NO FUNCIONA y se descarto:
  * la distancia minima al ideal mas cercano es un solo pico unimodal
    (moda 0.068) sin gap -> cualquier epsilon corta poblacion continua;
  * el ideal de Z15 cae casi en el centro de la nube de vidrio
    desordenado (el atomo mediano esta a 0.046 de Z15, contra 0.104 de
    hcp y 0.111 de Z14), asi que "el mas cercano gana" etiqueta ~70%
    del sistema como Z15 sin que esos atomos tengan nada de especial;
  * el ancho de la nube (sigma 3D ~ 0.077-0.104 segun nnn) es MAYOR
    que el epsilon propuesto.
La maquinaria de distancia (DIST_TYPES, distances_to_ideals) se deja
solo para que pilot_2e10_300K.py siga corriendo como evidencia de por
que se descarto. NO usarla para clasificar.

Estado por tipo despues de la revision:
  ico  -> criterio explicito del paper, validado en el piloto. SIN CAMBIOS.
  Z16  -> criterio del paper en dos etapas, con el corte de q4 sacado de
          los datos (ver mask_z16_q6window / mask_z16_revised).
  hcp,
  fcc  -> FUERA de la clasificacion dura. El paper tampoco los usa como
          curva numerica (Figs. 9-10 solo grafican n_ico y n_FK); quedan
          solo como puntos de referencia a dibujar sobre el mapa q4-q6
          (ver PLOT_ONLY_REFS).
  Z14,
  Z15  -> sin regla todavia. El paper no publica tolerancia para estos
          dos, asi que el criterio se decide DESPUES de mirar el mapa
          (q4, q6) por nnn (taskF_mapa_q4q6.py).
"""
import numpy as np

SPECIES = {1: "Zr", 2: "Cu", 3: "Ag"}
SPECIES_COLOR = {"Zr": "#2a78d6", "Cu": "#eb6834", "Ag": "#1baf7a"}

# indices 0-based de (Q4, Q6, W4hat, W6hat) dentro de cada fila del dump
NNN_COLS = {
    12: (2, 3, 4, 5),
    14: (6, 7, 8, 9),
    15: (10, 11, 12, 13),
    16: (14, 15, 16, 17),
}

# Tabla I, Ryltsev et al. 2016 (q4, q6, w4, w6) con el nnn propio de cada tipo
IDEAL = {
    "hcp": dict(nnn=12, q4=0.097,  q6=0.485,  w4=0.134,  w6=-0.012),
    "fcc": dict(nnn=12, q4=0.19,   q6=0.575,  w4=-0.159, w6=-0.013),
    "ico": dict(nnn=12, q4=1.4e-4, q6=0.663,  w4=-0.159, w6=-0.169),
    "Z14": dict(nnn=14, q4=0.015,  q6=0.406,  w4=0.134,  w6=-0.091),
    "Z15": dict(nnn=15, q4=0.039,  q6=0.316,  w4=0.134,  w6=-0.040),
    "Z16": dict(nnn=16, q4=0.0108, q6=0.1924, w4=0.1593, w6=0.013),
}

# hcp y fcc: SOLO para marcar sobre el mapa q4-q6 (nnn=12), no cuentan
# atomos. Decision de la Tarea F, alineada con el paper (Figs. 9-10 solo
# reportan n_ico y n_FK; hcp/fcc aparecen unicamente como referencia
# dibujada sobre el mapa).
PLOT_ONLY_REFS = ["hcp", "fcc"]

# DEPRECADO (ver nota de revision en el docstring): la clasificacion por
# bola de distancia se descarto. Se conserva solo para reproducir el
# diagnostico que la descarto.
DIST_TYPES = ["hcp", "fcc", "Z14", "Z15"]

# --- DECISION DE DISENO (no esta en el paper) -------------------------
# Ryltsev da umbrales explicitos SOLO para icosaedrico (q6 >= 0.6 y
# w6 <= -0.16, con nnn=12) y para Z16 (q4 <= 0.01 y |q6 - 0.1924| <=
# 0.01, con nnn=16). Para hcp, fcc, Z14 y Z15 no da numero, asi que
# aca se usa un criterio propio: distancia euclidea en el espacio
# (q4, q6, What6) al punto ideal de cada tipo, cada uno evaluado con
# SU PROPIO nnn (hcp/fcc con nnn=12, Z14 con nnn=14, Z15 con nnn=15).
# Se asigna el tipo mas cercano si esa distancia minima es < EPSILON;
# si ninguno cae dentro, el atomo queda "desordenado".
#
# Ojo con la sutileza: como cada tipo se mide con su propio nnn, las 4
# distancias NO viven en el mismo espacio de features. Es lo que pide
# la definicion (cada tipo ideal esta tabulado a su coordinacion), pero
# significa que "el mas cercano" compara distancias medidas sobre
# descriptores distintos del mismo atomo.
#
# W4hat queda FUERA de la distancia a proposito (el criterio es en
# (q4, q6, w6)), igual que en el criterio explicito del paper, que
# tampoco usa w4.
EPSILON_DEFAULT = 0.05
DIST_FEATURES = ("q4", "q6", "w6")


def load_boo_dump(path):
    """Devuelve el array (natoms, 18) de un *.boo.dump: id type + 16 cols."""
    with open(path) as fh:
        lines = fh.readlines()
    natoms = None
    start = None
    for i, line in enumerate(lines):
        if line.startswith("ITEM: NUMBER OF ATOMS"):
            natoms = int(lines[i + 1])
        elif line.startswith("ITEM: ATOMS"):
            start = i + 1
            break
    if natoms is None or start is None:
        raise ValueError(f"{path}: no encontre NUMBER OF ATOMS / ITEM: ATOMS")
    data = np.loadtxt(lines[start:start + natoms])
    if data.shape != (natoms, 18):
        raise ValueError(f"{path}: esperaba ({natoms}, 18), llego {data.shape}")
    return data


def features(data, nnn):
    """(q4, q6, w4hat, w6hat) como arrays, para un nnn dado."""
    c4, c6, cw4, cw6 = NNN_COLS[nnn]
    return data[:, c4], data[:, c6], data[:, cw4], data[:, cw6]


def distances_to_ideals(data, types=DIST_TYPES):
    """Distancia euclidea en (q4, q6, w6hat) de cada atomo al ideal de
    cada tipo, usando el nnn propio de ese tipo. Devuelve (natoms, ntypes)."""
    out = np.empty((data.shape[0], len(types)))
    for j, t in enumerate(types):
        ref = IDEAL[t]
        q4, q6, _w4, w6 = features(data, ref["nnn"])
        out[:, j] = np.sqrt(
            (q4 - ref["q4"]) ** 2 + (q6 - ref["q6"]) ** 2 + (w6 - ref["w6"]) ** 2
        )
    return out


def mask_ico(data):
    """Criterio EXPLICITO del paper: q6 >= 0.6 y w6 <= -0.16, con nnn=12."""
    _q4, q6, _w4, w6 = features(data, 12)
    return (q6 >= 0.6) & (w6 <= -0.16)


def mask_z16(data):
    """Criterio EXPLICITO del paper: q4 <= 0.01 y |q6 - 0.1924| <= 0.01,
    con nnn=16.

    OJO: en nuestras configs da CERO siempre. No es un bug: el minimo real
    de q4(nnn=16) sobre los 65536 atomos es 0.0129, por encima del ideal
    tabulado 0.0108, porque la distorsion termica a 300 K levanta q4 para
    todos. La mitad de q6 si la pasan ~2383 atomos. De ahi la version en
    dos etapas de abajo.
    """
    q4, q6, _w4, _w6 = features(data, 16)
    return (q4 <= 0.01) & (np.abs(q6 - 0.1924) <= 0.01)


def mask_z16_q6window(data, halfwidth=0.01):
    """Etapa 1 del Z16 revisado: la mitad de q6 del criterio del paper,
    |q6 - 0.1924| <= halfwidth con nnn=16. Esta mitad SI es alcanzable y
    es la que hace el trabajo de seleccion."""
    _q4, q6, _w4, _w6 = features(data, 16)
    return np.abs(q6 - IDEAL["Z16"]["q6"]) <= halfwidth


def mask_z16_revised(data, q4_cut, halfwidth=0.01):
    """Etapa 1 + etapa 2: dentro de la ventana de q6, exigir q4 <= q4_cut.

    q4_cut NO esta fijado aca a proposito: sale de mirar el histograma de
    q4 del subconjunto de la etapa 1 (taskF_z16_q4.py), no de un numero
    puesto a mano. El valor elegido y su justificacion se documentan al
    fijarlo.
    """
    q4, _q6, _w4, _w6 = features(data, 16)
    return mask_z16_q6window(data, halfwidth) & (q4 <= q4_cut)
