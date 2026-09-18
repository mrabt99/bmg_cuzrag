# Melt-quench Cu-Zr-Ag (N=65536) a 1e11 K/s — Negishi/SLURM

Segunda tasa de la comparación tipo Ryltsev del proyecto de tesis (Zr47.5 Cu47.5 Ag5 BMG).
La primera tasa (2e10 K/s) vive en el proyecto hermano `../65536_1e10/` — ver su
`README.md` y `ESTADO_PROYECTO_meltquench_1e10.md` para toda la historia previa: la
explosión de v1, la causa raíz (el potencial EAM no tiene core repulsivo válido por
debajo de ~1.2 Å), la decisión de no parchear el potencial, y el barrido de 7 tasas
(1e10, 2e10, 5e10, **1e11**, 2e11, 5e11, 1e12 K/s) sobre el tramo más caliente y
peligroso (2000→1900 K).

Este directorio es **autocontenido**: copia (no referencia por ruta relativa cruzada)
el potencial, el restart de partida a 2000 K, y el resultado del barrido a 1e11 K/s
desde `../65536_1e10/`, para que el proyecto no dependa de que ese otro directorio
siga existiendo en el mismo lugar.

## Por qué 1e11 K/s

En el barrido de tasas del proyecto anterior, 1e11 K/s **sobrevivió** el tramo más
caliente (2000→1900 K, 1 ns, 1000000 pasos) sin disparar el `fix halt` — de hecho con
más margen que 2e10 K/s (que fue la tasa más lenta que sobrevivió). Es la segunda tasa
natural para la comparación de propiedades (Tg, densidad, PE/átomo) en función de la
tasa de enfriamiento, con el mismo potencial sin tocar.

## Parámetros físicos (idénticos al proyecto 2e10)

- Zr47.5 Cu47.5 Ag5 BMG, N=65536 átomos, EAM/alloy `ZrCuAg.lammps.eam` (orden Zr Cu Ag),
  sin modificar.
- Ensemble NPT (Nose-Hoover), Tdamp=0.1 ps, **Pdamp=10.0 ps** ("v4").
- Timestep 1 fs, fijo (no bajar: rompería el análisis futuro de fonones/VDOS).
- `fix halt` de seguridad: T > 8000 K o |P| > ~10 GPa → corte duro.

## Los 17 segmentos, con su rango de T

Igual esquema que el proyecto 2e10: el segmento 0 (2000→1900 K) vive en `sweep_1e11/`
porque nació como parte del barrido de tasas; los 16 restantes (1900→300 K, 100 K
c/u) viven en `prod_1e11/`.

| # global | Rango T | Pasos (1e11 K/s) | Dónde vive | Restart de salida |
|---|---|---|---|---|
| 0 | 2000→1900 K | 1,000,000 (1 ns) | `sweep_1e11/` | `restart.sweep_1e11_1900K` |
| 1 | 1900→1800 K | 1,000,000 | `prod_1e11/` (seg01) | `restart.quench_N65536_1e11.seg1` |
| 2 | 1800→1700 K | 1,000,000 | `prod_1e11/` (seg02) | `restart.quench_N65536_1e11.seg2` |
| 3 | 1700→1600 K | 1,000,000 | `prod_1e11/` (seg03) | `...seg3` |
| 4 | 1600→1500 K | 1,000,000 | `prod_1e11/` (seg04) | `...seg4` |
| 5 | 1500→1400 K | 1,000,000 | `prod_1e11/` (seg05) | `...seg5` |
| 6 | 1400→1300 K | 1,000,000 | `prod_1e11/` (seg06) | `...seg6` |
| 7 | 1300→1200 K | 1,000,000 | `prod_1e11/` (seg07) | `...seg7` |
| 8 | 1200→1100 K | 1,000,000 | `prod_1e11/` (seg08) | `...seg8` |
| 9 | 1100→1000 K | 1,000,000 | `prod_1e11/` (seg09) | `...seg9` |
| 10 | 1000→900 K | 1,000,000 | `prod_1e11/` (seg10) | `...seg10` |
| 11 | 900→800 K | 1,000,000 | `prod_1e11/` (seg11) | `...seg11` |
| 12 | 800→700 K | 1,000,000 | `prod_1e11/` (seg12) | `...seg12` |
| 13 | 700→600 K | 1,000,000 | `prod_1e11/` (seg13) | `...seg13` |
| 14 | 600→500 K | 1,000,000 | `prod_1e11/` (seg14) | `...seg14` |
| 15 | 500→400 K | 1,000,000 | `prod_1e11/` (seg15) | `...seg15` |
| 16 | 400→300 K | 1,000,000 | `prod_1e11/` (seg16) | `...seg16` |

Aritmética: 100 K / 1e11 K/s = 1 ns = 1,000,000 pasos por tramo de 100 K (5x menos
pasos por tramo que a 2e10 K/s, consistente con ser 5x más rápida).

Relajación final: NPT 300 K / 0 Pa, 10 ns (1e7 pasos) en `prod_1e11/in.relax300K_1e11`
— misma duración/física que la relajación final del proyecto 2e10 (no depende de la
tasa que trajo al sistema hasta 300 K).

## Sobre el restart de partida — no es una "semilla"

`sweep_1e11/restart.sweep_1e11_1900K` (usado por `prod_1e11/seg01`) es el estado de
fase completo (posiciones + velocidades de los 65536 átomos + caja) al final del
tramo de barrido real que ya corrió en el proyecto anterior — **no** un número de
semilla RNG. Continuar la dinámica desde ahí es resumir la misma trayectoria física,
sin ningún `velocity create`; es determinista hacia adelante (misma versión de
LAMMPS/librerías). Mismo patrón que usó `prod_2e10/seg01` al continuar sin reseed
desde `restart.sweep_2e10_1900K`.

Un reseed de velocidades (`velocity create` con una semilla nueva) solo se agregaría
**reactivamente**, tramo por tramo, si algún segmento de esta cadena dispara el
`fix halt` — igual que se hizo en `prod_2e10/seg02` tras una explosión real — nunca
de antemano.

## El hueco del segmento 0: falta el dump de trayectoria

El script original del barrido (`sweep_1e11/in.ramp.lmp`, corrido en el proyecto
anterior como job 42877507) **no guardó dump de trayectoria** para el tramo
2000→1900 K — solo thermo (vía `fix ave/time`, en `sweep_1e11.log.dat`) y el restart
final. Mismo hueco que se identificó y tapó para `sweep_2e10/` en el proyecto
anterior (`in.ramp_with_dump.lmp` → `sweep_2e10_redump.log.dat`, verificado
determinista byte a byte contra el original).

`sweep_1e11/in.ramp_with_dump.lmp` repite exactamente la misma física (mismo restart
de entrada `../restart.equil2000K_N65536_1e10`, mismo potencial, mismo
Pdamp/Tdamp/dt, sin `velocity create`) agregando:
- `dump` cada 100,000 pasos (10 frames en el segmento) → `quench_N65536_1e11.dump.seg1_global`
- log de thermo con nombre distinto (`sweep_1e11_redump.log.dat`) para no pisar el original
- restart final con nombre distinto (`restart.sweep_1e11_1900K_redump`) — el restart
  que alimenta `prod_1e11/seg01` sigue siendo el original, no el redump

**Verificación pendiente tras correrlo:** comparar la última línea de
`sweep_1e11_redump.log.dat` contra `sweep_1e11.log.dat` (deben coincidir, o
prácticamente, salvo orden de suma en punto flotante) antes de dar el redump por
bueno — mismo chequeo que se hizo para `sweep_2e10_redump`.

## Estado actual / cómo correr

1. **`sweep_1e11/slurm_ramp_with_dump.sbatch`** — regenera el segmento 0 con dump
   (16 tasks, ~2h21m esperado según la corrida original, se pide 6h de margen).
   Independiente, sin dependencias — puede correr en paralelo con el resto.
2. **`prod_1e11/chain_submit.sh`** — encadena los 16 tramos de producción (1900→300K)
   + la relajación final, con `--dependency=afterok` (48 tasks, ~4h por tramo pedidas
   con margen sobre el ~1h esperado; 24h para el relax final de 10 ns).
   - `./chain_submit.sh` arranca desde seg01
   - `./chain_submit.sh --from N` retoma desde el tramo N (p.ej. tras una explosión)
3. **Recuperación ante explosión:** mismo procedimiento que en el proyecto 2e10 —
   cancelar el resto de la cadena, re-seedear velocidades en el último restart bueno
   con una semilla RNG nueva, agregar un hold NPT corto (0.5-2 ns) a esa T, retomar el
   tramo interrumpido, resubmitir el resto con `chain_submit.sh --from N`.

Ninguno de estos jobs se envió todavía a SLURM — quedan armados, pendientes de
confirmación antes de consumir horas de cómputo.

## Ver también

`../65536_1e10/README.md` y `../65536_1e10/ESTADO_PROYECTO_meltquench_1e10.md` para
la historia completa compartida (causa raíz del potencial, tabla completa del
barrido de 7 tasas, referencias Méndez & Jaroszewicz de Tg).
