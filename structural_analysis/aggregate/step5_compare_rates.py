#!/usr/bin/env python3
"""Tarea J - compare_rates. Un grafico por tipo de cluster, 3 tasas
superpuestas, un panel por especie central (equivalente Figs. 9-10 del
paper). n_ico_BOO y n_ico_Voronoi van en figuras SEPARADAS (no
comparables en escala, decision del usuario). Sobre abundancia_<tasa>.csv
(Tarea I), sin recalcular nada de BOO/Voronoi.

Mas la figura extra: razon n_ico(tasa mas lenta=2e10)/n_ico(tasa mas
rapida=1e12), en barras agrupadas por T, una figura para BOO y otra para
Voronoi (paneles por especie), para mostrar con la curva completa en T
la diferencia de sensibilidad que el piloto encontro solo en el punto de
300K.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RATES = ["2e10", "1e11", "1e12"]
RATE_COLOR = {"2e10": "tab:blue", "1e11": "tab:orange", "1e12": "tab:red"}
SPECIES = ["Zr", "Cu", "Ag"]
CLUSTER_COLS = ["n_ico_boo", "n_ico_voronoi", "n_Z14", "n_Z15", "n_Z16"]
CLUSTER_TITLE = {
    # ASCII simple a proposito (nada de Ŵ/What6 -- pedido del usuario,
    # para que el titulo se vea bien en cualquier fuente/render, tesis
    # incluida): W6_hat en vez del caracter unicode/mojibake.
    "n_ico_boo": "Icosaedrico (BOO: q6>=0.6, W6_hat<=-0.16, nnn=12)",
    "n_ico_voronoi": "Icosaedrico (Voronoi: indice <0,0,12,0>)",
    # CN total (Coordination de OVITO), NO indice de Voronoi exacto: no es
    # Frank-Kasper. Ver auditoria_Z_CN_vs_indice/INFORME.md.
    "n_Z14": "Coordinacion CN==14 (Voronoi, umbral de cara; no es Frank-Kasper)",
    "n_Z15": "Coordinacion CN==15 (Voronoi, umbral de cara; no es Frank-Kasper)",
    "n_Z16": "Coordinacion CN==16 (Voronoi, umbral de cara; no es Frank-Kasper)",
}

# piso minimo de conteo CRUDO (no fraccion) exigido en AMBAS tasas (2e10
# Y 1e12) antes de calcular una razon lenta/rapida -- sin esto, puntos
# con ~0-2 atomos en las dos tasas dan razones como 5.0 que no significan
# nada (ej. Cu a 1100K: BOO 2e10=10, 1e12=2 -> razon 5, puro ruido de
# conteo chico). Elegido mirando cuantos T (de 17) sobreviven por
# especie/metodo con distintos umbrales (10/15/20/25/30/50): 20 es el
# punto donde Zr-BOO y Ag-BOO (que nunca llegan a 20 en ninguna tasa)
# quedan en 0/17 sin ambiguedad, Cu-BOO se reduce a los 7 T mas frios
# (900K bajando) que es justo la zona con estadistica real, y
# Cu-Voronoi no pierde nada (17/17, nunca estuvo cerca del piso).
COUNT_FLOOR = 20

# ------------------------------------------------------------- carga
data = {}   # data[rate][especie][col] = (T_array, frac_array), T descendente
for rate in RATES:
    with open(f"abundancia_{rate}.csv") as f:
        rows = list(csv.DictReader(f))
    for sp in SPECIES:
        sub = [r for r in rows if r["especie"] == sp]
        sub.sort(key=lambda r: -int(r["T"]))
        T = np.array([int(r["T"]) for r in sub])
        for col in CLUSTER_COLS:
            data.setdefault(rate, {}).setdefault(sp, {})[col] = (
                T, np.array([float(r[col]) for r in sub]))

# conteos CRUDOS (no fraccion), solo para el piso de la figura de razon.
# vienen de tabla_clusters_51.csv, no de abundancia_<tasa>.csv (que ya
# esta en fraccion).
raw_counts = {}   # raw_counts[rate][especie][col] = (T_array, count_array)
RAW_COL = {"n_ico_boo": "N_ico_BOO", "n_ico_voronoi": "N_ico_Voronoi"}
with open("tabla_clusters_51.csv") as f:
    tabla_rows = list(csv.DictReader(f))
for rate in RATES:
    for sp in SPECIES:
        sub = [r for r in tabla_rows if r["tasa"] == rate and r["especie"] == sp]
        sub.sort(key=lambda r: -int(r["T_nominal"]))
        T = np.array([int(r["T_nominal"]) for r in sub])
        for col, rawcol in RAW_COL.items():
            raw_counts.setdefault(rate, {}).setdefault(sp, {})[col] = (
                T, np.array([int(r[rawcol]) for r in sub]))

# ------------------------------------------------ figuras por tipo de cluster
# OJO: hoy solo se regenera n_ico_boo. Las otras 4 (n_ico_voronoi, Z14,
# Z15, Z16) estan OBSOLETAS y sus PNG viven en obsoleto/: salen de
# abundancia_<tasa>.csv (pipeline A, step2), es decir Voronoi con
# relative_face_threshold=0.01 / edge_threshold=0.1 (los radios de Slater
# nunca se aplicaban: bug de compute(0)). Z14/Z15/Z16 son Coordination ==
# 14/15/16, NO el indice exacto <0,0,12,k>, asi que no son Frank-Kasper y
# no son comparables con tabla_fk_zr.csv. Version vigente (OVITO default,
# indice exacto): figs/compare_rates_*.png (step12). Ver
# auditoria_Z_CN_vs_indice/INFORME.md. Se excluyen del loop de escritura
# a proposito: no regenerarlas.
for col in ["n_ico_boo"]:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), sharex=True)
    for ax, sp in zip(axes, SPECIES):
        for rate in RATES:
            T, frac = data[rate][sp][col]
            ax.plot(T, frac * 100.0, "o-", ms=4, lw=1.8,
                   color=RATE_COLOR[rate], label=f"{rate} K/s")
        ax.set_title(sp)
        ax.set_xlabel("T (K)")
        ax.invert_xaxis()
        ax.grid(alpha=0.15)
    axes[0].set_ylabel("% de esa especie")
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle(CLUSTER_TITLE[col])
    fig.tight_layout()
    out = f"compare_{col}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"-> {out}")

# ------------------------------------------------ figura extra: razon lenta/rapida
SLOW, FAST = "2e10", "1e12"

def ratio_by_species_T(col):
    """Devuelve (T, razon, yerr) por especie. yerr viene de propagar el
    error de conteo tipo Poisson de cada N crudo (N_ico, no la fraccion):
    error relativo de la razon = sqrt(1/N_lenta + 1/N_rapida); yerr =
    razon * ese error relativo. Solo tiene sentido como aproximacion
    (asume Poisson puro, no correlacion entre configs ni entre especies),
    pero es la cota mas simple y estandar para "es un salto real o ruido
    de conteo" con lo que tenemos."""
    out = {}
    for sp in SPECIES:
        T_slow, f_slow = data[SLOW][sp][col]
        T_fast, f_fast = data[FAST][sp][col]
        _, n_slow = raw_counts[SLOW][sp][col]
        _, n_fast = raw_counts[FAST][sp][col]
        assert np.array_equal(T_slow, T_fast)
        enough = (n_slow >= COUNT_FLOOR) & (n_fast >= COUNT_FLOOR)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(enough, f_slow / np.where(f_fast > 0, f_fast, 1), np.nan)
            rel_err = np.sqrt(1.0 / np.where(n_slow > 0, n_slow, 1) +
                              1.0 / np.where(n_fast > 0, n_fast, 1))
        r = np.where(enough, r, np.nan)
        yerr = np.where(enough, r * rel_err, np.nan)
        out[sp] = (T_slow, r, yerr, n_slow, n_fast)
    return out

ico_boo_ratio = ratio_by_species_T("n_ico_boo")
ico_voro_ratio = ratio_by_species_T("n_ico_voronoi")

fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharex=True)
for ax, sp in zip(axes, SPECIES):
    T, r_boo, yerr_boo, _, _ = ico_boo_ratio[sp]
    _, r_voro, yerr_voro, _, _ = ico_voro_ratio[sp]
    width = 18
    valid_boo = ~np.isnan(r_boo)
    valid_voro = ~np.isnan(r_voro)
    ax.bar(T[valid_boo] - width/2, r_boo[valid_boo], width=width,
           yerr=yerr_boo[valid_boo], error_kw=dict(elinewidth=1.2, capsize=3,
           ecolor="#0b0b0b"), color="#4a3aa7", label="BOO")
    ax.bar(T[valid_voro] + width/2, r_voro[valid_voro], width=width,
           yerr=yerr_voro[valid_voro], error_kw=dict(elinewidth=1.2, capsize=3,
           ecolor="#0b0b0b"), color="#1baf7a", label="Voronoi")
    n_ok_boo = valid_boo.sum()
    n_ok_voro = valid_voro.sum()
    ax.set_title(sp)
    # nota de cobertura DENTRO del panel (no en el titulo, que a este
    # ancho de figura colisiona con el panel vecino): texto corto,
    # dos lineas, anclado en coordenadas de ejes para que no se escape.
    ax.text(0.03, 0.96, f"BOO: {n_ok_boo}/17 T\nVoronoi: {n_ok_voro}/17 T",
           transform=ax.transAxes, ha="left", va="top", fontsize=8,
           color="#52514e")
    ax.set_xlabel("T (K)")
    ax.invert_xaxis()
    ax.axhline(1.0, color="#898781", lw=1, ls=":")
    ax.grid(alpha=0.15, axis="y")
# OJO: esta figura esta OBSOLETA (parte Voronoi con umbral de cara); la
# version movida esta en obsoleto/. Reejecutar step5 la vuelve a escribir
# aqui: no usarla.
axes[0].set_ylabel(f"razon n_ico({SLOW}) / n_ico({FAST})")
axes[0].legend(frameon=False)
fig.suptitle(f"Sensibilidad a la tasa: razon icosaedrico {SLOW} (lenta) / "
             f"{FAST} (rapida), BOO vs Voronoi\n"
             f"(solo T con >= {COUNT_FLOOR} atomos en ambas tasas -- conteo "
             f"por panel dentro de cada grafico)", fontsize=11)
fig.tight_layout()
fig.savefig("compare_ico_ratio_slow_fast.png", dpi=140)
plt.close(fig)
print("-> compare_ico_ratio_slow_fast.png")

# resumen numerico en 300K (el punto que ya conociamos del piloto) + rango
print("\n=== razon 2e10/1e12 para Cu (especie dominante en ico), con error ===")
T, r_boo, yerr_boo, n_slow_boo, n_fast_boo = ico_boo_ratio["Cu"]
_, r_voro, yerr_voro, n_slow_voro, n_fast_voro = ico_voro_ratio["Cu"]
for i, t in enumerate(T):
    if np.isnan(r_boo[i]) and np.isnan(r_voro[i]):
        continue
    b = (f"BOO={r_boo[i]:.2f}+-{yerr_boo[i]:.2f} (N={n_slow_boo[i]}/{n_fast_boo[i]})"
         if not np.isnan(r_boo[i]) else "BOO=  --  ")
    v = (f"Voronoi={r_voro[i]:.2f}+-{yerr_voro[i]:.2f} (N={n_slow_voro[i]}/{n_fast_voro[i]})"
         if not np.isnan(r_voro[i]) else "Voronoi=  --  ")
    print(f"  T={t:5d}K  {b}   {v}")

# --- punto 3: el pico de Cu-BOO en 600K, se separa de 700K/500K? -------
print("\n=== punto 3: Cu-BOO, el pico de 600K se separa de sus vecinos? ===")
idx = {int(t): i for i, t in enumerate(T)}
for t_peak, t_left, t_right in [(600, 700, 500)]:
    i0, iL, iR = idx[t_peak], idx[t_left], idx[t_right]
    r0, e0 = r_boo[i0], yerr_boo[i0]
    rL, eL = r_boo[iL], yerr_boo[iL]
    rR, eR = r_boo[iR], yerr_boo[iR]
    print(f"  600K: {r0:.3f} +- {e0:.3f}  ->  intervalo [{r0-e0:.3f}, {r0+e0:.3f}]")
    print(f"  700K: {rL:.3f} +- {eL:.3f}  ->  intervalo [{rL-eL:.3f}, {rL+eL:.3f}]")
    print(f"  500K: {rR:.3f} +- {eR:.3f}  ->  intervalo [{rR-eR:.3f}, {rR+eR:.3f}]")
    overlap_L = (r0 - e0) <= (rL + eL)
    overlap_R = (r0 - e0) <= (rR + eR)
    print(f"  se superpone con 700K? {overlap_L}   se superpone con 500K? {overlap_R}")

# --- punto 4: cuantos T filtrados (N>=20) tienen barras que se solapan
# con AL MENOS UNO de sus vecinos inmediatos en T (mismo metodo/especie)?
print(f"\n=== punto 4: fraccion de T (ya filtrados por N>={COUNT_FLOOR}) cuya "
      f"barra de error se solapa con algun T vecino ===")
for label, ratio_dict in [("BOO", ico_boo_ratio), ("Voronoi", ico_voro_ratio)]:
    for sp in SPECIES:
        Tsp, rsp, esp, _, _ = ratio_dict[sp]
        valid_idx = [i for i in range(len(Tsp)) if not np.isnan(rsp[i])]
        if len(valid_idx) < 2:
            print(f"  {label:8s} {sp}: solo {len(valid_idx)} T valido(s), "
                  f"no hay vecino para comparar")
            continue
        n_overlap = 0
        for k, i in enumerate(valid_idx):
            lo_i, hi_i = rsp[i] - esp[i], rsp[i] + esp[i]
            neighbors = []
            if k > 0:
                neighbors.append(valid_idx[k - 1])
            if k < len(valid_idx) - 1:
                neighbors.append(valid_idx[k + 1])
            overlaps = any(
                (lo_i <= rsp[j] + esp[j]) and (hi_i >= rsp[j] - esp[j])
                for j in neighbors
            )
            n_overlap += overlaps
        print(f"  {label:8s} {sp}: {n_overlap}/{len(valid_idx)} T se solapan "
              f"con al menos un vecino inmediato")
