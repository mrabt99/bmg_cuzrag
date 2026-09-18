#!/usr/bin/env python3
"""Tarea L - agregado final por regimen para compare_ico_ratio_slow_fast.png.

CORRECCION DE METODO respecto de lo que quedo insinuado en la Tarea K:
NO promediar las 5 razones punto a punto del regimen vitreo. Se suman
los CONTEOS CRUDOS N_ico de todos los T de cada regimen (por separado
para 2e10 y 1e12), y se calcula UNA razon sobre esos conteos agrupados:

    razon_regimen = N_ico_pool(2e10) / N_ico_pool(1e12)
    error = razon_regimen * sqrt(1/N_pool_2e10 + 1/N_pool_1e12)

(la formula de error es la misma sqrt(1/N1+1/N2) de la Tarea K, aplicada
a los conteos YA sumados, no combinando los errores de cada T por
separado -- promediar errores punto a punto habria sido el error de
metodo que el usuario corrigio).

Definicion de regimen: Tg propio de cada tasa (promedio de Tg_V y Tg_Ep
de Tg_resultados.csv, Etapa 1) +-100K:
    liquido/liq. subenfriado: T > Tg + 100
    vitreo:                   T < Tg - 100
    (la banda intermedia, ni liquido ni vitreo, se excluye del agregado)

Tg_2e10 = (646.1+637.2)/2 = 641.65 K  -> liquido: 1900-800K (12 T),
                                          vitreo: 500-300K (3 T)
Tg_1e12 = (703.2+697.9)/2 = 700.55 K  -> liquido: 1900-900K (11 T),
                                          vitreo: 600-300K (4 T)

CORRECCION (Tarea M): la version con ventana "propia" (Tg+-100K de CADA
tasa por separado) se descarto como entregable, para las 4 combinaciones
que faltaban (Zr y Ag, Voronoi) igual que ya se habia hecho para Cu en
la Tarea L. No es un problema exclusivo de Cu: como Tg_2e10 != Tg_1e12,
la ventana "propia" compara T ABSOLUTAS distintas entre tasas (ej.
liquido de 2e10 llega a 800K, el de 1e12 se corta en 900K) -- eso rompe
la comparacion de fondo, porque Ryltsev compara ambas tasas a la MISMA
T, no a la misma distancia relativa a su propio Tg. La version "propia"
se deja solo como diagnostico impreso (para que quede constancia de
por que se descarto), pero YA NO se escribe en tabla_regimen_pool.csv.
El entregable final usa SIEMPRE la ventana COMPARTIDA (interseccion de
T entre 2e10 y 1e12: liquido=1900-900K, vitreo=500-300K), igual para
las 4 combinaciones (Cu-BOO, Zr-Voronoi, Cu-Voronoi, Ag-Voronoi).

Advertencia para el texto de la tesis (Tarea L, punto 4): el error
sqrt(1/N1+1/N2) asume conteos de Poisson INDEPENDIENTES entre atomos.
El orden icosaedrico es espacialmente correlacionado -- Ryltsev et al.
describen "superclusters" de icosaedros conectados a tasas lentas, que
por construccion NO son eventos independientes (un icosaedro extra
suele aparecer pegado a otro, no en un lugar al azar). Este error es
por lo tanto un PISO (subestima la incertidumbre real), no un intervalo
de confianza riguroso.
"""
import csv
import math

SPECIES = ["Zr", "Cu", "Ag"]
RATES = ["2e10", "1e12"]
BUFFER = 100.0

TG = {  # promedio Tg_V/Tg_Ep, Tg_resultados.csv (Etapa 1)
    "2e10": (646.1 + 637.2) / 2,
    "1e12": (703.2 + 697.9) / 2,
}

with open("tabla_clusters_51.csv") as f:
    rows = list(csv.DictReader(f))


def regime_Ts(rate, regime):
    Tg = TG[rate]
    Ts = sorted({int(r["T_nominal"]) for r in rows if r["tasa"] == rate},
                reverse=True)
    if regime == "liquido":
        return [t for t in Ts if t > Tg + BUFFER]
    elif regime == "vitreo":
        return [t for t in Ts if t < Tg - BUFFER]
    raise ValueError(regime)


print("=== paso 1: T por regimen y tasa (buffer +-100K sobre Tg propio) ===")
for rate in RATES:
    print(f"  {rate} (Tg={TG[rate]:.1f}K): liquido={regime_Ts(rate,'liquido')} "
          f"({len(regime_Ts(rate,'liquido'))} T)  "
          f"vitreo={regime_Ts(rate,'vitreo')} ({len(regime_Ts(rate,'vitreo'))} T)")

COL = {"BOO": "N_ico_BOO", "Voronoi": "N_ico_Voronoi"}


def pool(rate, especie, metodo, Ts):
    col = COL[metodo]
    return sum(int(r[col]) for r in rows
               if r["tasa"] == rate and r["especie"] == especie
               and int(r["T_nominal"]) in Ts)


def ratio_err(n_slow, n_fast):
    if n_slow == 0 or n_fast == 0:
        return float("nan"), float("nan")
    r = n_slow / n_fast
    err = r * math.sqrt(1.0 / n_slow + 1.0 / n_fast)
    return r, err


targets = [("Cu", "BOO")] + [(sp, "Voronoi") for sp in SPECIES]

header = f"{'especie':8s} {'metodo':8s} {'regimen':9s} {'N_pool_2e10':>12s} {'N_pool_1e12':>12s} {'razon':>8s} {'error':>8s}"

print("\n=== diagnostico DESCARTADO: ventana propia de cada tasa (no se guarda) ===")
print(header)
for especie, metodo in targets:
    for regime in ["liquido", "vitreo"]:
        T_2e10 = regime_Ts("2e10", regime)
        T_1e12 = regime_Ts("1e12", regime)
        n_slow = pool("2e10", especie, metodo, T_2e10)
        n_fast = pool("1e12", especie, metodo, T_1e12)
        r, err = ratio_err(n_slow, n_fast)
        print(f"{especie:8s} {metodo:8s} {regime:9s} {n_slow:12d} {n_fast:12d} "
              f"{r:8.3f} {err:8.3f}")

# --- entregable final: ventana COMPARTIDA (interseccion de ambas), para
# las 4 combinaciones (Cu-BOO + Zr/Cu/Ag-Voronoi) ----------------------
print("\n=== paso 2/3: tabla final -- ventana COMPARTIDA (interseccion 2e10/1e12) ===")
print(header)
final_rows = []
for especie, metodo in targets:
    for regime in ["liquido", "vitreo"]:
        T_shared = sorted(set(regime_Ts("2e10", regime)) &
                          set(regime_Ts("1e12", regime)), reverse=True)
        n_slow = pool("2e10", especie, metodo, T_shared)
        n_fast = pool("1e12", especie, metodo, T_shared)
        r, err = ratio_err(n_slow, n_fast)
        print(f"{especie:8s} {metodo:8s} {regime:9s} {n_slow:12d} {n_fast:12d} "
              f"{r:8.3f} {err:8.3f}   (T compartidos: {T_shared})")
        final_rows.append({"especie": especie, "metodo": metodo, "regimen": regime,
                           "N_pool_2e10": n_slow, "N_pool_1e12": n_fast,
                           "razon": round(r, 4), "error": round(err, 4)})

with open("tabla_regimen_pool.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(final_rows[0].keys()))
    w.writeheader()
    w.writerows(final_rows)
print(f"\n-> tabla_regimen_pool.csv ({len(final_rows)} filas, solo ventana compartida)")
