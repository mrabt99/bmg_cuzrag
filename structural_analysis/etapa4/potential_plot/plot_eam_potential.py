#!/usr/bin/env python3
"""Grafica el potencial EAM/alloy usado en todo el proyecto (ZrCuAg.lammps.eam, formato
setfl DYNAMO, generado por Howard Sheng/GMU). Parsea F(rho), rho(r) y r*phi(r) por par
directamente del archivo (sin libreria eam especial, solo texto plano).

Entrada : /scratch/negishi/aaabt/bmg_cuzrag/rate_2e10/ZrCuAg.lammps.eam
Salida  : fig_eam_potential.png (3 paneles: F(rho), rho(r), phi(r))
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

POT_FILE = Path("/scratch/negishi/aaabt/bmg_cuzrag/rate_2e10/ZrCuAg.lammps.eam")
OUT = Path(__file__).resolve().parent.parent / "fig_eam_potential.png"

lines = POT_FILE.read_text().splitlines()
# 3 lineas de comentario, luego "Nelements el1 el2 el3", luego "Nrho drho Nr dr cutoff"
elements = lines[3].split()[1:]
nelem = len(elements)
nrho, drho, nr, dr, rcut = lines[4].split()[:5]
nrho, nr = int(nrho), int(nr)
drho, dr, rcut = float(drho), float(dr), float(rcut)
print(f"elementos: {elements}  nrho={nrho} drho={drho}  nr={nr} dr={dr}  rcut={rcut}")

# el resto son todos numeros, con 1 linea de cabecera "Z mass a0 lattice" antes de cada
# bloque F_i / rho_i. Concatenamos todos los floats de las lineas restantes y navegamos por
# indice, saltando esas lineas de cabecera por separado (no son floats puros, tienen la
# palabra del lattice al final -> las manejamos aparte).
idx = 5
F = {}
RHO = {}
for el in elements:
    header = lines[idx].split()
    idx += 1
    # el resto del archivo desde aca es todo numeros: juntamos lineas hasta tener lo necesario
    nvals_needed = nrho + nr
    vals = []
    while len(vals) < nvals_needed:
        vals.extend(float(x) for x in lines[idx].split())
        idx += 1
    F[el] = np.array(vals[:nrho])
    RHO[el] = np.array(vals[nrho:nrho + nr])

PHI = {}
pair_order = []
for i in range(nelem):
    for j in range(i + 1):
        pair_order.append((elements[i], elements[j]))
for (ei, ej) in pair_order:
    vals = []
    while len(vals) < nr:
        vals.extend(float(x) for x in lines[idx].split())
        idx += 1
    r_phi = np.array(vals[:nr])
    PHI[(ei, ej)] = r_phi  # tabulado como r*phi(r), unidades eV*Angstrom

rho_grid = np.arange(nrho) * drho
r_grid = np.arange(nr) * dr
r_grid_safe = np.where(r_grid == 0, np.nan, r_grid)

colors = {"Zr": "#1f77b4", "Cu": "#d62728", "Ag": "#7f7f7f"}
pair_colors = plt.get_cmap("tab10")

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

axF, axRho, axPhi = axes
for el in elements:
    axF.plot(rho_grid, F[el], color=colors.get(el), lw=1.6, label=el)
axF.set_xlabel(r"$\rho$")
axF.set_ylabel(r"$F(\rho)$, eV")
axF.set_title("Funcion de embedding")
axF.legend(frameon=False)
axF.axhline(0, color="0.8", lw=0.7, zorder=0)

for el in elements:
    axRho.plot(r_grid, RHO[el], color=colors.get(el), lw=1.6, label=el)
axRho.set_xlabel(r"$r$, $\mathrm{\AA}$")
axRho.set_ylabel(r"$\rho(r)$")
axRho.set_title("Densidad electronica atomica")
axRho.set_xlim(0, rcut)
axRho.legend(frameon=False)

for k, (ei, ej) in enumerate(pair_order):
    phi = PHI[(ei, ej)] / r_grid_safe  # deshace el r*phi(r) tabulado
    label = f"{ei}-{ej}"
    axPhi.plot(r_grid, phi, color=pair_colors(k), lw=1.6, label=label)
axPhi.set_xlabel(r"$r$, $\mathrm{\AA}$")
axPhi.set_ylabel(r"$\phi(r)$, eV")
axPhi.set_title("Potencial de pares")
axPhi.set_xlim(1.5, rcut)
axPhi.set_ylim(-1.0, 2.0)
axPhi.axhline(0, color="0.8", lw=0.7, zorder=0)
axPhi.legend(frameon=False, fontsize=8)

fig.suptitle(f"Potencial EAM/alloy -- {'-'.join(elements)} (Sheng/GMU, {POT_FILE.name})", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUT, dpi=200)
print(f"escrita: {OUT}")
