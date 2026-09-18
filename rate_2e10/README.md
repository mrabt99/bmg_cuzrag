# Melt-quench Cu-Zr-Ag (N=65536) a 2e10 K/s — Negishi/SLURM

Repo de backup completo: scripts, logs de LAMMPS/SLURM, curvas termodinámicas `.log.dat`,
este README, `ESTADO_PROYECTO_meltquench_1e10.md`, y los binarios (restarts + dumps de
trayectoria de la cadena de producción). Ver "Qué se puede y qué no se puede reproducir con
esto" al final para el detalle de cobertura y limitaciones reales de los datos.

Ver **ESTADO_PROYECTO_meltquench_1e10.md** para la historia completa del proyecto (v1 explotó,
v2 descartado, v3 en sol67, causa raíz del potencial, decisión de tasa, sweep, cadena final).
Este README es la referencia operativa: qué corrió, con qué parámetros, y qué archivo es cada
cosa.

## TL;DR

- Zr47.5 Cu47.5 Ag5 BMG, N=65536 átomos, EAM/alloy `ZrCuAg.lammps.eam` (orden Zr Cu Ag),
  **sin modificar** (ver por qué abajo).
- El potencial **no tiene core repulsivo válido por debajo de ~1.2 Å** (4 de los 6 pares son
  atractivos ahí) — causa raíz de las explosiones estocásticas de v1 y del sweep. No depende
  de la tasa de enfriamiento ni del barostato; depende de cuánto tiempo pasa el sistema
  líquido y caliente. Detalle y datos en `check_eam_core.py` + `ESTADO_PROYECTO...md`.
- Se decidió **no parchear el potencial** (invalidaría la comparación entre tasas tipo
  Ryltsev) y en cambio medir empíricamente la tasa más lenta que sobrevive el tramo más
  caliente (2000→1900 K).
- Resultado del barrido (2000→1900 K, un solo seed por tasa, mismo restart de partida):

  | Tasa | Duración del tramo | Resultado | Explotó a (ns abs.) | % del tramo |
  |---|---|---|---|---|
  | 1e10 K/s | 10 ns | explotó | 1.188 ns | 11.9% |
  | 2e10 K/s | 5 ns | **sobrevivió entero** | — | — |
  | 5e10 K/s | 2 ns | explotó | 1.014 ns | 50.7% |
  | 1e11 K/s | 1 ns | sobrevivió | — | — |
  | 2e11 K/s | 0.5 ns | sobrevivió | — | — |
  | 5e11 K/s | 0.2 ns | sobrevivió | — | — |
  | 1e12 K/s | 0.1 ns | sobrevivió (control) | — | — |

  Patrón no monótono (1e10 falla, 2e10 sobrevive, 5e10 falla) — evidencia de que es un evento
  estocástico ligado al tiempo absoluto a alta T, no un umbral determinístico de tasa. Toda la
  corrida final es a **2e10 K/s constante**, 2000→300 K, potencial original sin tocar.
- Referencia Méndez & Jaroszewicz para este sistema: **Tg≈692 K** (criterio energía potencial),
  **Tg≈664 K** (Wendt-Abraham modificado), σy≈2.3 GPa, UTS≈2.7 GPa.
- Ensemble: **NPT** (Nose-Hoover), Tdamp=0.1 ps, **Pdamp=10.0 ps** ("v4" — no 1.0 como la
  referencia a 1e12 K/s del profesor).
- Timestep 1 fs, fijo (no bajar a 0.5 fs: rompería el análisis futuro de fonones/VDOS por
  autocorrelación de velocidades).
- **17 segmentos de 100 K** cubren 2000→300 K: segmento 1 = `sweep_2e10/` (2000→1900 K, 5 ns,
  5e6 pasos), segmentos 2-17 = `prod_2e10/` (1900→300 K, 5 ns c/u). **No son 16, son 17** —
  el segmento 1 vive en `sweep_2e10/` porque nació como parte del barrido de tasa, no de
  `prod_2e10/`.
- Relajación final: NPT 300 K / 0 Pa, 10 ns, en `prod_2e10/in.relax300K_2e10` (aparte del
  templado, mismo motivo que el profesor: agregarla sin reiniciar LAMMPS dispara error con
  `fix ave/time` ya activo — ver gotchas).
- Solo hizo falta **un reintento** en toda la cadena de producción: el segmento 2 (1800→1700 K)
  explotó la primera vez y se recuperó con semilla de velocidades nueva + 0.5 ns de hold NPT
  antes de retomar la rampa (ver `prod_2e10/in.quench_2e10_seg02.lmp`, que documenta el
  reintento en sus propios comentarios). Los segmentos 3-17 corrieron limpios a la primera.

## Los 17 segmentos, con su rango de T

| # global | Rango T | Dónde vive | Archivo de restart |
|---|---|---|---|
| 1 | 2000→1900 K | `sweep_2e10/` | `restart.sweep_2e10_1900K` |
| 2 | 1900→1800 K | `prod_2e10/` (seg01) | `restart.quench_N65536_2e10.seg1` |
| 3 | 1800→1700 K | `prod_2e10/` (seg02, reintento con nueva semilla + hold) | `restart.quench_N65536_2e10.seg2` |
| 4 | 1700→1600 K | `prod_2e10/` (seg03) | `restart.quench_N65536_2e10.seg3` |
| 5 | 1600→1500 K | `prod_2e10/` (seg04) | `restart.quench_N65536_2e10.seg4` |
| 6 | 1500→1400 K | `prod_2e10/` (seg05) | `restart.quench_N65536_2e10.seg5` |
| 7 | 1400→1300 K | `prod_2e10/` (seg06) | `restart.quench_N65536_2e10.seg6` |
| 8 | 1300→1200 K | `prod_2e10/` (seg07) | `restart.quench_N65536_2e10.seg7` |
| 9 | 1200→1100 K | `prod_2e10/` (seg08) | `restart.quench_N65536_2e10.seg8` |
| 10 | 1100→1000 K | `prod_2e10/` (seg09) | `restart.quench_N65536_2e10.seg9` |
| 11 | 1000→900 K | `prod_2e10/` (seg10) | `restart.quench_N65536_2e10.seg10` |
| 12 | 900→800 K | `prod_2e10/` (seg11) | `restart.quench_N65536_2e10.seg11` |
| 13 | 800→700 K | `prod_2e10/` (seg12) | `restart.quench_N65536_2e10.seg12` |
| 14 | 700→600 K | `prod_2e10/` (seg13) | `restart.quench_N65536_2e10.seg13` |
| 15 | 600→500 K | `prod_2e10/` (seg14) | `restart.quench_N65536_2e10.seg14` |
| 16 | 500→400 K | `prod_2e10/` (seg15) | `restart.quench_N65536_2e10.seg15` |
| 17 | 400→300 K | `prod_2e10/` (seg16) | `restart.quench_N65536_2e10.seg16` |

Segmentos 14-15 (700→500 K) son los que cruzan la banda de Tg de referencia (664-692 K) —
zona más sensible para chequear cristalización o cualquier anomalía estructural.

## Estructura de carpetas

```
.
├── README.md, ESTADO_PROYECTO_meltquench_1e10.md   <- documentación
├── ZrCuAg.lammps.eam                                <- potencial (texto, setfl), sin modificar
├── check_eam_core.py                                <- diagnóstico del core repulsivo faltante
├── check_pressure.py                                <- monitor de T/P sobre los .log.dat*
├── chain_submit.sh                                   <- encadena jobs SLURM (afterok)
│
├── sweep_1e10/, sweep_2e10/, sweep_5e10/,            <- barrido de tasa, 2000->1900K,
│   sweep_1e11/, sweep_2e11/, sweep_5e11/, sweep_1e12/    un job independiente por tasa
│   (sweep_2e10/ = SEGMENTO 1 de la cadena final, el único que sobrevivió y se continuó)
│
├── prod_2e10/                                        <- SEGMENTOS 2-17: producción real,
│   in.quench_2e10_seg01.lmp ... seg16.lmp                 1900->300K a 2e10 K/s constante
│   in.relax300K_2e10                                      (la numeración interna de estos
│   slurm_seg01.sbatch ... seg16.sbatch, slurm_relax.sbatch  archivos es 1-16, pero en la
│                                                            cadena global son los segmentos
│                                                            2-17: seg01 acá = segmento 2 global)
│
├── in.quench_1e10*, in.relax300K_1e10, seg01-17.lmp  <- INTENTO ABANDONADO a 1e10 K/s puro
│   (explotó, ver ESTADO_PROYECTO "v1"/"RESULTADO DEL BARRIDO"). Se deja como referencia
│   histórica de por qué se pasó a 2e10, no se usó para la muestra final.
│
└── restart.melt_N65536                               <- punto de partida (2000K, de in.melt
                                                            del profesor), referenciado por
                                                            sweep_*/ y por el intento a 1e10
```

## Por qué Pdamp=10 y no 1.0

Un intento previo (v1, en sol67) usó Pdamp=1.0 (igual que la referencia a 1e12 K/s del
profesor) y **explotó a 1.67 ns** por lo que en ese momento se interpretó como runaway del
barostato. Subir Pdamp a 10.0 ps le da al barostato más inercia sin dejar de permitir que la
caja contraiga libremente a P≈0 (a diferencia de NVT, que evita el síntoma pero deja tensión
interna sin relajar). **Corrección posterior (05/09):** el diagnóstico más profundo mostró que
el barostato es víctima, no causa — el defecto real está en el core del potencial. Pdamp=10 se
mantuvo de todos modos en toda la corrida final (no hace daño, y da algo más de margen), pero
el verdadero control de la explosión terminó siendo la elección de tasa (2e10), no Pdamp.

## LAMMPS: versión y gotchas de sintaxis encontrados acá

- La corrida de referencia del profesor (1e12 K/s) se hizo con LAMMPS (12 Dec 2018). En
  Negishi solo hay módulos `lammps/20220623` y `lammps/2022.5.04`. Se usó **lammps/20220623**
  (`module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623`). No es reproducción bit-a-bit.
- **Gotcha 1:** un `print` con variables equal-style (`v_Tinst`, `v_PEat`, etc.) **después** de
  `write_restart`/`write_data`/`write_dump` tira `ERROR: Compute temp thermo_temp used in
  variable thermo keyword between runs is not current` en esta versión (no pasaba en 12 Dec
  2018). Fix: esos `print` van **antes** de `write_restart`.
- **Gotcha 2:** `reset_timestep 0` con un `fix ave/time` ya activo tira `ERROR: Cannot reset
  timestep with time-dependent fix ave/time defined` (variante del mismo problema que
  `reset_timestep` + barostato NPT reportado en sol67, mensaje distinto según versión). Fix:
  sacar el `reset_timestep` cosmético, no hace falta para la física.
- Potencial: el input original del profesor apuntaba a `../ZrCuAg.lammps.eam`, que en su
  carpeta real NO existe un nivel arriba. Acá se copió el potencial junto con los scripts.

Cuenta SLURM usada en todos los `.sbatch`: **`strachan`**. Layout de recursos: 1 nodo, 48 tasks
en la cadena de producción (`prod_2e10/`), 16 tasks en el barrido (`sweep_*/`, por cupo
compartido del grupo — ver `AssocGrpGRES` en el historial de la sesión si hace falta reproducir
esa parte).

## Cómo reproducir / continuar

1. `sweep_2e10/in.ramp.lmp` + `sweep_2e10/slurm_ramp.sbatch`: segmento 1 (2000→1900 K).
   Ya corrido, restart final: `sweep_2e10/restart.sweep_2e10_1900K` (incluido en este repo).
2. `prod_2e10/chain_submit.sh` encadena `in.quench_2e10_seg01.lmp` (arranca del restart de
   sweep_2e10) hasta `seg16.lmp` + `in.relax300K_2e10`, con `sbatch --dependency=afterok`
   (nunca `afterany`: si un segmento muere, el siguiente no debe leer un restart malo o
   inexistente). Soporta `--from N` para retomar después de una falla.
3. Si un segmento explota: cancelar la cola, editar ese `in.quench_2e10_segNN.lmp` para
   agregar `velocity all create <T> <semilla_nueva> mom yes rot yes dist gaussian` después del
   `read_restart` + un hold NPT breve (0.3-0.5 ns) a esa T antes de retomar la rampa, correr un
   test de sintaxis corto, relanzar con `chain_submit.sh --from N`. Ejemplo real de esto:
   `prod_2e10/in.quench_2e10_seg02.lmp`.
4. Monitorear con `check_pressure.py <log.dat...>` (detecta columnas por nombre exacto, ojo
   que `title2` no etiqueta la columna de paso que `fix ave/time` siempre antepone — hay que
   sumarla a mano al indexar si se reusa este script en otro lado).

## Cuidados que siguen aplicando

- Vigilar caída abrupta e inusual de energía potencial durante el enfriamiento (más allá de la
  caída suave esperada en la transición vítrea) — señal de cristalización espontánea (Ryltsev
  et al. 2016 reportan cristalitos de Cu2Zr a 1.5e9 K/s).
- **Esto es un solo seed.** El protocolo está validado (17 segmentos limpios con un solo
  reintento) pero la estadística no — confirmar con 2-3 seeds más antes de tratar esta muestra
  como completamente caracterizada.

## Qué se puede y qué no se puede reproducir con esto

| Cubierto por este backup | NO cubierto |
|---|---|
| Protocolo completo rearmable: los 17 inputs + sbatch de cada segmento, potencial verificado (checksum), causa raíz del defecto documentada (`check_eam_core.py` + `ESTADO_PROYECTO...md`) | VDOS por autocorrelación de velocidades sobre el enfriamiento: todos los dumps (producción y el segmento 1 regenerado) muestrean cada 500.000 pasos (500 ps entre frames) — **demasiado grueso para autocorrelación de velocidades**. Hace falta una corrida NVE dedicada con salida fina (cada pocos fs) arrancando desde un restart puntual; esa corrida no existe todavía en ningún lado |
| Trayectoria completa de los 17 segmentos, incluido el segmento 1 (2000→1900 K): el sweep original no guardaba dump ahí, se regeneró aparte (`sweep_2e10/in.ramp_with_dump.lmp`, mismo restart de partida, sin re-samplear nada) y se verificó que el estado final coincide byte a byte con el sweep original — `sweep_2e10/quench_N65536_2e10.dump.seg1_global` | Difusión por MSD durante el enfriamiento — mismo problema de resolución temporal que el VDOS |
| Curva V/T y Tg vía criterio de volumen o energía potencial — 17 puntos de 100 K con thermo completo (`*.log.dat*`) | Comparación estadística entre semillas: esto es una sola trayectoria; la tasa 2e10 no está confirmada con réplicas |
| BOO (q4/q6/w4/w6), RDF, clasificación de clusters/Voronoi en cada uno de los 17 puntos de 100 K — los restarts dan la configuración exacta (posiciones+velocidades) en cada uno | |
| Reproducir el protocolo completo tramo por tramo, incluida la recuperación real del segmento 3 (nueva semilla + hold, documentada en los comentarios de `prod_2e10/in.quench_2e10_seg02.lmp`) | |

Restarts y dumps de la cadena de producción (`prod_2e10/`) SÍ están en este repo (no en Drive
ni en un `.gitignore` — esa versión anterior del plan se descartó).

`quench_N65536_1e10.dump.validate_seed2` (200MB, la trayectoria de la validación que corrió
limpia 4.5 ns con la semilla buena en el intento abandonado a 1e10 — no hace falta para
reproducir la metodología final a 2e10, queda como evidencia de esa corrida) superaba el
límite de 100MB/archivo de GitHub sin Git LFS, así que está partido en 3 con `split -b 90M`:
`quench_N65536_1e10.dump.validate_seed2.part_aa/ab/ac`. Para reconstruirlo:

```
cat quench_N65536_1e10.dump.validate_seed2.part_* > quench_N65536_1e10.dump.validate_seed2
```
