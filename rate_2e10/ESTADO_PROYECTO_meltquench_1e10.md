# Estado del proyecto — Melt-quench Cu-Zr-Ag a 1e10 K/s (sol67)

Handoff traído por Agustín desde su sesión de Claude Code en sol67 (03/09/2026). Documenta el estado real de la corrida de enfriamiento a 1e10 K/s, incluida la "explosión de Abt". Reemplaza en la práctica al prompt genérico `claude/protocolo-quench-1e10.md`, que fue escrito sin este contexto y contiene consejos equivocados (ver nota al final).

## Objetivo de la corrida

Zr47.5Cu47.5Ag5, N=65536 átomos, potencial EAM/alloy `ZrCuAg.lammps.eam` (orden Zr-Cu-Ag). Comparar la estructura del vidrio a 1×10^10 K/s vs. la referencia 1×10^12 K/s. Referencia Méndez & Jaroszewicz: Tg≈692 K (criterio PE), 664 K (Wendt-Abraham modificado), σy≈2.3 GPa, UTS≈2.7 GPa.

## Cluster

Ojo: la tesis BMG normalmente corre en **sol67 / sol125** (LabSim, CNEA, scheduler **PBS**), pero para ESTA corrida de 1e10 Agustín quiere usar **Negishi** (Purdue, scheduler **SLURM**) por única vez. El profe le dejó la referencia 1e12 en `/scratch/negishi/aaabt/65536`. Así que el nuevo input/job de 1e10 se arma para Negishi/SLURM (`sbatch`, `--dependency=afterok`), NO con la receta PBS/qsub de sol67. Las lecciones de física (Pdamp=10, tensión del NVT, ~100× de pasos) valen igual en cualquier cluster; lo que sigue de sol67 es solo referencia histórica de dónde ocurrió la saga v1-v4.

### Referencia histórica sol67 (donde ocurrió la explosión)

sol67/sol125, scheduler PBS (OpenPBS/OpenHPC), NO SLURM pese a la wiki de LabSim.

- Cola activa: `long` (workq deshabilitada). Comandos `qsub`/`qstat`/`pbsnodes -a`, solo desde head node.
- LAMMPS: `module load lammps/2025` (OpenMPI 4.1.4 + Intel).
- Producción: clase c128, N=64 cores físicos (~5.7× más rápido en wall-time absoluto que c40/N=16). Nombres de clase = threads, no cores físicos (c8→4, c40→20, c48→24, c128→64); `--bind-to core` falla más allá del límite físico.
- Selección de clase va en la línea de qsub: `-l select=1:ncpus=64:c128=True`, no dentro del `.pbs`.
- Receta paralela funcional: `module load lammps/2025` → `export OMP_NUM_THREADS=1` → hostfile manual (`sort -u $PBS_NODEFILE | awk '{print $1" slots=N"}' > hostfile`) → `mpiexec --mca ras ^tm --bind-to core --map-by core --report-bindings -hostfile hostfile -n N lmp -in file.in -log log.txt`. Los cuatro flags son obligatorios.

## Historia de versiones

- **v1 `~/meltquench_1e10/` — EXPLOTÓ (01/09/2026, "explosión de Abt").** NPT a 1e10 K/s. Runaway catastrófico del barostato a 1.67 ns en el líquido a 2000 K. Causa raíz confirmada: 10 ns/segmento a 1e10 da ~100× más pasos que la referencia 1e12 (100 ps/segmento), lo que dispara la probabilidad de una fluctuación térmica espontánea → fuerza repulsiva EAM masiva → feedback positivo del barostato NPT. Densidad/Lx planos hasta el instante del colapso. Restart mal equilibrado y resonancia del barostato descartados definitivamente. No se escribieron restarts antes de la explosión.
- **v2 `~/meltquench_1e10_v2/` — descartado.** Bajó dt a 0.5 fs (manteniendo NPT). Descartado por Jaroszewicz: comprometería el análisis futuro de fonones/VDOS vía autocorrelación de velocidades, que necesita dt=1 fs.
- **v3 `~/meltquench_1e10_v3/` — activo (NVT, lanzado 02/09).** Cadena de ~19 jobs en c128/64 cores: equilibración NVT 2000K/2ns → 17 rampas NVT 100K/10ns c/u (afterok) → NVT breve + NVE 10ns final. Todos con `fix halt` (T>8000K, `error hard`) y tracking de presión `v_Pinst` en `fix log2` como criterio de validación acordado.
  - **Riesgo de v3:** a volumen fijo el sistema no puede contraerse durante el enfriamiento. La contracción volumétrica real es 9.22%. Tensión estimada a 300 K: ~7-13 GPa (B~80-120 GPa × 9.22%), 3-5× σy/UTS. Riesgo real de fractura/cavitación NO cubierto por `fix halt` (que solo dispara por temperatura).
  - Monitoreo: `~/check_pressure.py` parsea los `.log.dat`, detecta columnas por nombre de header (fix del falso positivo `P(` dentro de `Ep(`), convierte atm→GPa, compara contra `--fracture` (default 2.3 GPa). Primer check real (paso 3.04M tramo 1, T~1970K): P=-0.150 GPa, ~6.5% del margen. Ventana crítica: tramos 13-14 (~800-600 K, cerca de Tg).
- **v4 (NPT, Pdamp=10.0) — propuesto por Jaroszewicz (02/09 noche), no arrancado.** Volver a NPT pero subir Pdamp de 1.0 a 10.0 ps (solo el barostato; Tdamp queda en 0.1 ps). Idea: 10× más inercia del barostato promedia la presión en una ventana más larga e ignora el pico repulsivo transitorio del EAM que disparó la explosión, rompiendo el feedback positivo. Ventaja sobre v3: el sistema contrae libremente a P≈0, elimina el problema de tensión por volumen fijo. Cuidado: la explosión apareció a 1.67 ns, así que un smoke test de 2000 pasos NO alcanza; requiere una corrida de validación de ~4-5 ns (2-3× el tiempo a la explosión) antes de comprometer la cadena completa. Pendiente: construir `in.quench_1e10` con Pdamp=10 partiendo del NPT original (no de v3).
- **Control de ensamble `~/control_ensamble_1e12/` — diseñado, no lanzado.** `in.control_NVT_1e12`: NVT a 1e12 que sigue la curva real V(T) del NPT 1e12 vía `change_box remap` en cada uno de los 17 bordes de segmento. Aísla el efecto puro de ensamble (NPT vs. NVT) sobre Tg. Fix requerido: `run 0 post no` antes del primer `change_box`.

## CAUSA RAÍZ ENCONTRADA (05/09/2026): el potencial no tiene core repulsivo

Diagnóstico hecho sobre `ZrCuAg.lammps.eam` (setfl, nr=1001, dr=0.0065, cutoff=6.5). Descartados antes: listas de vecinos (`Dangerous builds = 0`, `neighbor 2.0 bin`, `neigh_modify delay 0 every 1 check yes`) y la función de embedding (F(rho) monótona creciente pasando el mínimo en los tres elementos: Zr min en rho=0.628, Cu en 0.725, Ag en 0.588).

**El problema está en φ(r) por debajo de ~1.2 Å. Cuatro de los seis pares son ATRACTIVOS ahí:**

| par | φ(0.4 Å) | φ(0.8 Å) | φ(1.0 Å) | φ/ZBL a 0.8 Å |
|---|---|---|---|---|
| Zr-Zr | +687.5 | +276.6 | +153.3 | 0.56 |
| Cu-Zr | **−296.2** | **−97.4** | **−42.0** | −0.25 |
| Cu-Cu | +53.3 | +29.4 | +20.6 | 0.09 |
| Ag-Zr | **−207.4** | **−62.5** | **−22.7** | −0.11 |
| Ag-Cu | **−240.8** | **−78.4** | **−33.3** | −0.18 |
| Ag-Ag | **−41.5** | +17.3 | +27.9 | 0.03 |

φ creciente con r significa fuerza atractiva. Un par Cu-Zr que entra adentro de ~1.2 Å se atrae cada vez más fuerte al acercarse: es un pozo sin fondo. Los dos pares que sí son repulsivos son inservibles igual por lo blandos (Cu-Cu da 2% de ZBL a 0.4 Å).

**Esto explica todas las observaciones:** el colapso de un par libera energía sin cota (de ahí 5 MeV y 638.832 K); el evento es raro y aleatorio porque depende de que una fluctuación empuje un par adentro del radio crítico (tiempos dispersos 0.71/1.19/1.67 ns, dependencia con la semilla); no hay precursor macroscópico (todo plano hasta el instante); Pdamp no cambia nada porque el barostato es la víctima y no la causa; y la referencia a 1e12 nunca explotó solo porque corría 100 ps por tramo en vez de 10 ns, o sea 100× menos exposición al mismo defecto latente.

**Conclusión:** no es un problema de tasa de enfriamiento, ni de ensamble, ni del termostato. El potencial nunca tuvo core repulsivo válido y recién se ve al correr lo suficientemente largo. Las muestras existentes a 1e12 no quedan invalidadas (no hubo colapso), pero cualquier corrida larga necesita el core arreglado.

Dos arreglos posibles: (a) empalmar core ZBL a las funciones de par del setfl (blend de Fermi centrado en rc=1.6 Å, b=12 Å⁻¹, transición ~1.1-2.1 Å, muy por debajo del primer vecino ~2.5-2.8 Å; script `patch_eam_zbl.py`); (b) enfriar rápido (1e12) el tramo caliente y lento (1e10) solo la región relevante. Ver decisión abajo.

## Decisión (05/09/2026): tasa variable, NO parchear el potencial

Agustín (ya hablado con Jaroszewicz) decidió NO parchear el potencial, para no arriesgar la comparación entre tasas del estudio tipo Ryltsev (todas las tasas deben usar el mismo potencial). En su lugar: **enfriar rápido (1e12 K/s) el tramo caliente 2000→1500 K, y lento (1e10 K/s) solo de 1500→300 K.**

Fundamento físico: la zona peligrosa del potencial es la de alta T (más energía térmica → más cruces de la barrera hacia el pozo atractivo), y es justo donde la tasa de enfriamiento NO afecta al vidrio final, porque el líquido arriba de ~1500 K está en equilibrio y se re-equilibra rápido. La estructura del vidrio la define la región del líquido subenfriado y la transición vítrea (Tg~680 K), que se recorre a 1e10. Así la muestra es legítimamente "1e10 K/s".

Implementación (la más barata): arrancar del restart de la referencia 1e12 a 1500 K (ese tramo ya está enfriado a 1e12, no se re-corre) y seguir 1500→300 K a 1e10 en cadena de 12 jobs afterok, NPT Pdamp=10, potencial original. Alternativa: correr el propio tramo 2000→1500 a 1e12 desde restart.equil2000K.

Riesgo abierto y cómo se mitiga: NO está verificado que los tramos intermedios/bajos (por debajo de 1500 K) sean seguros a 1e10, porque ninguna corrida sobrevivió más allá del tramo 1 para saberlo. Por eso ANTES de la cadena larga se hace una prueba de estabilidad: correr ~2 ns a 1e10 desde los restarts de la referencia 1e12 a ~1200 K y ~800 K. Si no explotan, el peligro está confinado al tramo caliente y el plan es sólido. Si explotan, hay que replantear el punto de cruce.

Refinamiento (05/09, más tarde): antes de fijar el cruce en 1500 K, se decidió MEDIR primero. Barrido en paralelo (jobs independientes, no encadenados) del primer tramo 2000→1900 K (el más peligroso), todos desde el mismo restart seed2, variando la tasa (1e10, 2e10, 5e10, 1e11, 2e11, 5e11, 1e12) para encontrar la más lenta que no explota. Según el resultado se usará esa tasa constante (más limpio para la comparación, sin punto de cruce que justificar) o una tasa variable con cruce. Ojo metodológico: el evento es estocástico, así que un solo job por tasa da una muestra; la tasa elegida se confirma después con 2-3 seeds más.

**RESULTADO DEL BARRIDO 2000→1900 K (06/09/2026), corrido en Negishi (`/scratch/negishi/aaabt/65536_1e10`):**

| Tasa | Duración rampa | Resultado | Explotó a (ns abs) | % de rampa |
|---|---|---|---|---|
| 1e10 K/s | 10 ns | explotó | 1.188 | 11.9% |
| 2e10 K/s | 5 ns | **SOBREVIVIÓ** (tramo completo) | — | — |
| 5e10 K/s | 2 ns | explotó | 1.014 | 50.7% |
| 1e11 K/s | 1 ns | sobrevivió | — | — |
| 2e11 K/s | 0.5 ns | sobrevivió | — | — |
| 5e11 K/s | 0.2 ns | sobrevivió | — | — |
| 1e12 K/s | 0.1 ns | sobrevivió (control) | — | — |

Lo revelador: las dos que fallaron explotaron a ~1.0-1.2 ns absolutos (T aún ~1950-1990 K), casi el mismo instante pese a ser 12% de la rampa para 1e10 y 51% para 5e10. Confirma que el hazard está ligado al tiempo absoluto a alta T, NO a la velocidad de la rampa: proceso estocástico con escala característica ~1-1.5 ns, concentrado en alta T. **2e10 es la tasa más lenta viable para el tramo caliente** (1e10, la única más lenta probada, falló). Job sweep_2e10: NPT Pdamp=10, 16 cores, dt=1 fs, 5M pasos (3e6→8e6), llegó a 1900 K limpio. Este job (2000→1900 K) es el **segmento 1** de la cadena completa a 2e10.

## Decisión (06/09/2026): continuar a 2e10 CONSTANTE hasta 300 K

Como 2e10 sobrevivió el tramo MÁS caliente y los tramos más fríos son estrictamente más seguros (hazard cae con T), se descarta el cruce de tasa y el barrido por ventana. Plan: continuar `sweep_2e10` directo desde su restart final (1900 K, ~paso 8e6) a **2e10 constante 1900→300 K**, en 16 segmentos adicionales de 100 K (5M pasos = 5 ns c/u), NPT Pdamp=10, mismo potencial y dt=1 fs, restart cada 100 K y dump cada 10 K (como la referencia 1e12), cadena afterok. Sumado al segmento 1 (`sweep_2e10`, 2000→1900 K), la cadena completa 2000→300 K queda en **17 segmentos**, no 16 (16 es solo la parte de producción `prod_2e10`, 1900→300 K). La muestra queda legítimamente "2e10 K/s", limpia para la comparación tipo Ryltsev sin punto de cruce que justificar.

Riesgo residual: los segmentos aún templados (1900→~1400 K) acumulan algo de peligro; se mitiga con restarts frecuentes y **validando primero SOLO el segmento 1900→1800** antes de encadenar el resto. Si ese segmento explota, es la señal para volver a tasa variable en la zona caliente. Pendiente aún: confirmar la tasa elegida con 2-3 seeds más, y vigilar caída abrupta de PE (cristalización Ryltsev, cristalitos Cu2Zr) durante toda la cadena. Prompt de continuación para el asistente en Negishi entregado (06/09).

## CADENA 2e10 COMPLETA (11/09/2026): terminó sin explotar

La cadena completa de **17 segmentos** (segmento 1 = `sweep_2e10`, 2000→1900 K; segmentos 2-17 = `prod_2e10`, 1900→300 K de 100 K en 100 K) terminó. El gráfico de volumen/átomo vs. T sobre el set completo, 2000→300 K, muestra una curva limpia y monótona, con inflexión visible justo en la banda 664-692 K donde caen las Tg de referencia (PE 692 K, Wendt-Abraham 664 K) — consistente con lo esperado, sin caídas anómalas de volumen que delaten cristalización. Primera muestra completa de melt-quench Cu-Zr-Ag a tasa lenta (2e10 K/s) que sobrevive entera con el potencial sin parchear.

Solo hizo falta un reintento en toda la cadena de producción: el segmento 2 (1800→1700 K) explotó en su primer intento (paso 15938000, 2.94 ns dentro del tramo, ~1750 K) y se recuperó con semilla de velocidades nueva + hold NPT breve de 0.5 ns a esa T antes de retomar la rampa. Los segmentos 3-17 corrieron limpios a la primera. El job de relajación final (300 K/10 ns) falló al instante por un gotcha de sintaxis (`reset_timestep 0` con `fix ave/time` ya activo — ver gotchas) y se relanzó sin ese `reset_timestep` cosmético.

Pendiente real: esto es un solo seed. La confirmación con 2-3 seeds adicionales (mencionada arriba) sigue sin hacerse; el resultado de hoy valida el protocolo, no cierra la estadística.

Siguiente paso (11/09/2026): backup de los archivos necesarios para reproducir la metodología tipo Ryltsev, separando GitHub (scripts/logs/README, texto) de Google Drive (restarts y datos binarios). Prompt de backup entregado al asistente en Negishi, con README que debe listar los 17 segmentos (no 16): segmento 1 del barrido (2000→1900) + 16 de producción (1900→300).

## Gotchas de LAMMPS registrados

- `reset_timestep 0` + `run ... pre no` → `ERROR: Virial was not tallied on needed timestep` con barostato NPT. Sacar el `reset_timestep` cosmético.
- `reset_timestep 0` con `fix ave/time` ya activo (sin `pre no`, en LAMMPS 20220623/Negishi) → `ERROR: Cannot reset timestep with time-dependent fix ave/time defined`. Mismo gotcha en el fondo (reset_timestep chocando con un fix con estado temporal ya activo), mensaje distinto según versión. Mismo arreglo: sacar el `reset_timestep` cosmético, no hace falta para la física.
- `change_box` justo después de `read_restart` → `ERROR: Change_box parameter not allowed after reading restart file with per-atom info`. Fix: `run 0 post no` primero.
- `fix ave/time`: `P(` es substring de `Ep(`, falso positivo en detección de columnas. Detectar por nombre exacto.
- Cadenas donde cada job lee el restart del anterior: usar `depend=afterok` (no afterany).

## Nota sobre el prompt viejo

`claude/protocolo-quench-1e10.md` se escribió sin este contexto. Corrección importante: dice "Tdamp/Pdamp no dependen de la tasa, no tocar" cuando Pdamp=10 es justamente el fix de v4, así que ese consejo es al revés de lo que necesita. En cambio la parte de Negishi/SLURM (`sbatch`, `--dependency`) SÍ aplica: esta corrida va a Negishi por decisión de Agustín. Usar este doc como fuente de verdad para el resto.
