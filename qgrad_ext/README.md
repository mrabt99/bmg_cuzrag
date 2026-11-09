# qgrad_ext — rama de qgrad que extiende la tasa lenta por debajo de Tg (plan A)

**Estado: LANZADA el 2026-09-30 (ver "Jobs").** El usuario eligió el plan A y pidió partir seg10 en 10a/10 (restart en 550 K sin costo extra de cómputo).
Ubicación: `bmg_cuzrag/qgrad_ext/` (al lado de `rate_gradiente_gammamin/`). No se modificó nada de `qgrad/`, `rate_gradiente_gammamin/`, las cadenas 2e10/1e11/1e12 ni `structural_analysis/`
(solo se LEE de ahí).

## Motivo
qgrad (2e10→1.5e9 K/s entre 1500 y 700 K, luego 1e11 K/s entre 700 y 300 K) quedó con Ep(300 K) solo −0.65 meV/at por debajo de 2e10, aunque Cu-ico BOO da 1.19 (5σ). Hipótesis: el cambio a
tasa rápida en 700 K ocurre por encima de la Tg de la tasa lenta (sistema todavía líquido), y el tramo rápido "rejuvenece" la energía. La rama sigue lenta hasta bien por debajo de esa Tg.

## Paso 0 — auditoría de qgrad (lo que la rama copia EXACTO)
| | valor (idéntico en los 12 segmentos de qgrad) |
|---|---|
| dt | 1 fs (`timestep 0.001`, "fijo, no tocar") |
| ensemble | NPT `iso 0.0 0.0`, Tdamp 0.1 ps, Pdamp 10.0 ps |
| potencial | `eam/alloy` `ZrCuAg.lammps.eam` (Zr Cu Ag), `neighbor 2.0 bin` (copia con sha256 verificado) |
| fix halt | cada 1000 pasos: `v_Tinst > 8000` y `v_Pinst > 1.0e5`, ambos `error hard`. **Nota:** es `Pinst > 1e5`, no `abs(...)` (el comentario dice |P|): solo dispara con presión positiva |
| frecuencias tramo lento (seg08) | thermo 66700, `ave/time 667 100 66700` (~1000 filas de log.dat), dump cada 6 666 667 pasos (11 frames con velocidades) |
| frecuencias tramo rápido (seg09–12) | thermo 1000, `ave/time 10 100 1000`, dump cada 100 000 |
| restart | solo al final de cada segmento; sin `velocity create` (cada tramo continúa del anterior) |
| tasas | seg08 (800→700 K) = **1.500e9 K/s** (66 666 667 pasos); seg09–12 (700→300 K) = **1.000e11 K/s** (1 000 000 pasos c/u) |

Restart de 700 K (`restart.quench_N65536_gradmin.seg08` → copia `restart.start_700K_from_qgrad_seg8`, sha256 `e31cdd9b…`): verificado contra el log —el paso 0 de seg09 coincide en todos los
dígitos con el "terminado" de seg08: T=700.7221 K, P=−504.765 bar, Lx=105.3409 Å, ρ=7.346249 g/cc, Ep=−4.778764 eV/at.
Velocidad medida en qgrad: seg08 289 pasos/s (48 ranks); tramos rápidos 272–281; relajación 10 ns = 8 h 47 min.

## Paso 1 — dónde cortar la tasa lenta (`analisis/tg_extrapolacion.py`, `analisis/tg_extrapolacion.csv`)
Tg de Etapa 1 (media de Tg(V) y Tg(Ep), ventanas estándar): 2e10 → 641.7 K, 1e11 → 666.1 K, 1e12 → 700.5 K; pendiente 34.6 K/década (casi perfectamente lineal en log γ).
**Extrapolación a 1.5e9 K/s: Tg = 602.8 K.** Solo 3 puntos; la extrapolación (1.125 décadas por debajo de la tasa más lenta) **asume linealidad en log γ**. El error estándar del ajuste (±0.4 K, 1 grado de libertad)
NO es la incertidumbre real: la real es la sensibilidad a las ventanas de ajuste, **577–664 K** (media de Tg(V) y Tg(Ep); 564–668 K separando). La Tg medida del perfil qgrad es 631.9 K (perfil mixto).
**T_switch = Tg − 100 K = 502.8 K → 500 K (elegido, plan A).** Opción −50 K → 550 K (no elegida; ahorraría ~31.5 h).

## Paso 2 — cadena (plan A). Segmentos, misma numeración que qgrad (segNN termina a la misma T que en qgrad)
| tramo | T ini→fin | tasa (K/s) | pasos | horas est. | walltime pedido | restart de salida |
|---|---|---|---|---|---|---|
| seg09 | 700→600 | 1.5e9 | 66 666 667 | 64.1 | 5-16:19:00 | `restart.quench_N65536_gradext.seg09` |
| seg10a | 600→550 | 1.5e9 | 33 333 333 | 32.1 | 2-19:54:00 | `…seg10a` (restart en 550 K: punto de corte opcional) |
| seg10 | 550→500 (=T_switch) | 1.5e9 | 33 333 334 | 32.1 | 2-19:54:00 | `…seg10` |
| seg11 | 500→400 | 1e11 | 1 000 000 | 1.0 | 02:33:00 | `…seg11` |
| seg12 | 400→300 | 1e11 | 1 000 000 | 1.0 | 02:33:00 | `…seg12` |
| relax | 300 K, 10 ns NPT | — | 10 000 000 | 8.8 | 1-00:00:00 | `restart.glass300K_N65536_gradext` |
En 10a y 10 se conservan las MISMAS frecuencias de thermo/ave-time/dump que seg08 de qgrad (misma ventana de promedio de 66.7 ps): ~500 filas de log.dat y 6 frames por tramo de 33.3 M pasos.
Total **145.3 M pasos ≈ 139 h ≈ 5.8 días** de cómputo (128.2 h lentos + 2.0 h rápidos + 8.8 h relajación), a 289 pasos/s en los lentos. (Referencia: la cadena de qgrad tardó 8.0 días de calendario para 191.8 h de
cómputo, es decir sin esperas apreciables.) 48 ranks, `-A strachan -p cpu`, `--mpi=pmi2`. Walltime = 2.2× la estimación a 300 pasos/s (mismo criterio que qgrad).
Las configuraciones de 1900 a 700 K (01–13 de la grilla) son las mismas que las de qgrad; la rama aporta 14–17 (600, 500, 400, 300 K).
Riesgos: la cuenta tiene un tope de grupo `gres/hp_cpu=320` (un job de 48 cores puede quedar en cola); un fallo en seg09/seg10 pierde hasta ~64 h.
Recuperación ante explosión: la de qgrad, SIN hold isotérmico (reseed de velocidades en el restart de entrada del tramo; si falla >2–3 veces, parar y reportar).

**Archivos:** `in.quench_ext_seg09/10a/10/11/12.lmp`, `in.relax300K_ext` (generados por `build_inputs.py` a partir de los de qgrad; solo cambian rampa de T, nombres de archivo `rate=gradext` y restart de entrada —ver `diff`),
`slurm_seg09/10a/10/11/12.sbatch`, `slurm_relax.sbatch`, `chain_submit.sh`, `restart.start_700K_from_qgrad_seg8`, `ZrCuAg.lammps.eam`, `monitor_ext.py`.

## Paso 3 — monitoreo (`monitor_ext.py`, sin análisis estructural)
Por tramo: borde exacto (terminado del tramo previo vs paso 0 del siguiente) y de ventanas (tolerancias de la auditoría de qgrad: ΔT ≤ 3 K, ΔLx ≤ 0.013 Å, Δρ ≤ 0.002, ΔEp ≤ 0.3 meV/at),
stderr limpio (sin ruido PMIx/munge), sin halt/ERROR/lost atoms, `check_crystallization.py`; figura `monitor/monitor_Ep_gradext.png` (Ep(T), Ep−Ep(2e10), pendiente local) con 2e10, 1e11, 1e12 y qgrad.
**Alarma** (solo informa; **la decisión de cortar es del usuario**): pendiente local dEp/dT en 2 ventanas consecutivas de 50 K dentro de ±5 % de 0.143 meV/at/K (pendiente vítrea de las 4 curvas en 300–500 K;
la líquida es 0.255) con deriva de la brecha Ep−Ep(2e10) < 0.3 meV/at. Probado sobre los logs de qgrad (todo OK; brecha a 315 K = −0.66 meV/at). Con seg10 partido en 10a/10 hay restart en 550 K (fin de seg10a): si la alarma
dispara, se puede cortar ahí y lanzar el tramo rápido desde `restart.quench_N65536_gradext.seg10a` (decisión del usuario).

## Predicciones (escritas ANTES de correr; asumen Ep y ln N lineales en log10 γ; pendiente Ep 4.23 meV/at por década; errores Poisson piso, ±0.04 en Zr-FK)
| escenario | Ep(300 K, relajado) − Ep(2e10) | Zr-FK Z14 / Z15 / Z16 vs 2e10 | Cu-ico BOO / Voronoi vs 2e10 |
|---|---|---|---|
| 1: el tramo rápido era lo que borraba (tasa efectiva ≈ 1.5e9) | ≈ −4.8 meV/at (≥ 3) | ≈ 1.07 / 1.14 / 1.19 (empiezan a separarse de 2e10) | la ley predice ≈ 1.5 / 1.2 |
| 2: los icosaedros de Cu ya sobrevivían | ≈ −3 a −5 meV/at | tienden a 1.07–1.19 | ≈ 1.19 / 1.10 (poco cambio respecto de qgrad; cambia sobre todo la energía) |
| 3: nada cambia respecto de qgrad | ≈ −0.65 ± 0.5 meV/at | ≈ 1.0 | ≈ 1.19 / 1.10 → la ventana 1500–700 K no alcanzó; el límite no es el tramo rápido |
Referencias (qgrad vs 2e10, vidrio 500/400/300 K): Ep −0.65 meV/at; Zr Z14/Z15/Z16 = 1.033/1.018/1.015; Cu-ico BOO 1.193, Voronoi 1.105. (Comparación futura con el pipeline de `structural_analysis/`
copiando el criterio de `qgrad/`: manifest con las configs 14–17 de esta rama; las 01–13 ya existen.)

## Jobs
Lanzados el 2026-09-30 con `./chain_submit.sh` (encadenados con `afterok`; `squeue -u $USER`):
| tramo | job |
|---|---|
| seg09 (700→600 K, 1.5e9) | 43612277 |
| seg10a (600→550 K, 1.5e9) | 43612278 |
| seg10 (550→500 K, 1.5e9) | 43612279 |
| seg11 (500→400 K, 1e11) | 43612280 |
| seg12 (400→300 K, 1e11) | 43612281 |
| relax (300 K, 10 ns) | 43612282 |
Si un tramo falla la cadena se corta sola: relanzar con `./chain_submit.sh --from <tramo>`.

## Resultado y pipeline estructural (2026-10-06)
Cadena completada (jobs 43612277-82, todos COMPLETED). `monitor_ext.py` (con `--start-log` = log de seg08 de qgrad): bordes exactos OK, sin halt/cristalización; el único "REVISAR" (seg11, borde de ventanas dLx −0.019 Å) es ruido de una ventana de 1 ps (σ_Lx 0.008 Å).
Relajación 300 K: Ep −4.84084 eV, ρ 7.4770, deriva −0.001 meV/at/ns. **Ep − Ep(2e10) = −3.66 meV/at** (qgrad: −0.64).
Estructura (`manifest_ext.txt`: 14=600 K, 15=500, 16=400, 17r=300 relajado, 17u=300 sin relajar; jobs ext_dump 43852734, ext_boo 43852735, ext_voro 43852759, ext_post 43852780; scripts en `aggregate/`: `voronoi_ext.py`, `consolidate_ext.py`, `regimen_ext.py`, `fig_ext.py`; salidas `tabla_compare_rates_5curvas.csv`, `tabla_regimen_ext.csv`, `tabla_porT_ext.csv`, `figs/compare_ext_5curvas.png`).
Razón vítrea (500/400/300 K, ventana fija igual a la de qgrad) ext/2e10: Cu-ico BOO 1.263±0.040, Voronoi 1.171±0.017, Zr Z14/Z15/Z16 = 1.041/1.159/1.066. ext/qgrad: Cu-ico BOO 1.059, Voronoi 1.060, Z15 1.138 (4.6σ), Z14 1.007, Z16 1.049.
Lectura: el tramo lento hasta 500 K sube Ep (−3.0 meV extra) y Zr-Z15 (≈ +14 % sobre qgrad), pero Cu-ico solo +6 % sobre qgrad y Z14/Z16 poco; la predicción del escenario 1 (Z14/15/16 ≈ 1.07/1.14/1.19, Cu-ico ≈ 1.5/1.2 vs 2e10) se cumple solo en parte (Z15 sí, Z16 y Cu-ico BOO no llegan a la ley log-lineal). Una sola muestra por tasa; errores = piso Poisson.
**No corrido:** RDF/BADF/Cu-Cu, Laves, Tg del perfil ext, validación del 5º color (el validador de la skill dataviz no corre con el node del cluster).

## Cu-Cu y Laves (2026-10-07)
**Cu-Cu 2º pico a 300 K, 10 frames de la relajación** (`etapa4/cucu_multiframe_ext.py`, copia del de qgrad con qgrad_ext agregado; reproduce exacto los números de las 4 curvas; salidas `etapa4/cucu_multiframe_ext.csv`, `…resumen_ext.csv`, `fig_R_CuCu_multiframe_5curvas.png`; job 43883256).
| | profundidad media de frames (± sem) | del g promediado (± bootstrap) |
|---|---|---|
| 1e12 | 0.321 ± 0.004 | 0.285 ± 0.007 |
| 1e11 | 0.353 ± 0.003 | 0.318 ± 0.005 |
| 2e10 | 0.369 ± 0.005 | 0.329 ± 0.006 |
| qgrad | 0.376 ± 0.004 | 0.344 ± 0.007 |
| **qgrad_ext** | **0.387 ± 0.004** | **0.350 ± 0.006** |
qgrad_ext − 2e10 = +0.0185 ± 0.0066 (2.8σ) / +0.021 ± 0.008 (2.6σ); qgrad_ext − qgrad ≈ +0.012 ± 0.006 (≈ 2σ). Orden monótono en la lentitud, la rama queda en el extremo. Errores = SOLO ruido térmico entre frames (una muestra por tasa). OJO: la validación contra LAMMPS solo vale para 1e12 (único cuyo último frame = restart); para los demás Δg ≤ 0.25 es esperable (frame ≠ restart final).
**Laves** (`laves/`, `laves_analysis_ext.py` = copia con el set qgrad_ext, criterios IDÉNTICOS, features/ = symlinks de solo lectura; job 43883270): **diamante en subred Zr = 0** en todas las configs; Zr Laves-like (T&B) = 0 en 06–17r de la rama (qgrad 17u tenía 1); Z16 con 4 Zr vecinos = 9 (17r). Sin evidencia de Laves; cluster máximo Laves-like ≤ 1. Salidas `laves/laves_conteos_ext.csv`, `laves_clusters_ext.csv`, `fig_laves_conteos_5curvas.png`.
**Sigue sin hacerse:** RDF/BADF completos de las configs 14–17 (solo Cu-Cu 300 K), Zr-Cu cruzado, clusters Z16+ico, Tg del perfil de la rama.

## Comparación con la detección de Ryltsev (2026-10-07, `ryltsev_cmp/`, job 43883326)
Ryltsev et al. no usan criterio cuantitativo: (i) abundancias BOO de ico/Z16 con salto en γmin, (ii) mapa Fig. 11a (centros ico coloreados por q6) donde ven "un área con orden regular", (iii) "análisis detallado de snapshots" → Cu2Mg (Z16 de Zr + ico de Cu con caras compartidas). Sistema: Cu64.5Zr35.5, N=5000, Mendelev EAM, 1800→300 K.
Hecho acá (300 K relajado, `ryltsev_cmp.py`, salidas `ryltsev_cmp.csv`, `fig11a_style_300K.png` = loncha de 12 Å con todos los átomos por q6):
| | N ico (BOO) | factor de agrupación de ico (vecinos < 4.5 Å / azar) | Z16 BOO estricto | Z16 BOO solo ventana q6 | cluster mayor Zr-Z16+Cu-ico (Voronoi) |
|---|---|---|---|---|---|
| 2e10 | 770 | 4.90 | 0 | 2383 | 30 |
| 1e11 | 575 | 5.00 | 0 | 2284 | 29 |
| qgrad | 903 | 4.73 | 0 | 2340 | 51 |
| qgrad_ext | 952 | 4.60 | 0 | 2267 | 49 |
Lectura: la rama tiene más centros ico pero NO más agrupados que 2e10 (el factor baja con la lentitud); cluster mayor ≈ qgrad (49 vs 51, ~1.4 % de esos poliedros, antes ya descartado como grano). El Z16 BOO estricto del paper da 0 a 300 K (distorsión térmica; ver classify/boo_types.py), así que la comparación de Z16 usa Voronoi. En el mapa se ven acumulaciones de ico en qgrad y qgrad_ext, sin patrón de red a simple vista; una loncha 2D no puede mostrar cristalinidad. Limitaciones: el factor de agrupación incluye la correlación de empaquetamiento y no distingue red periódica de cúmulos amorfos; no hay semejanza alguna con su grano medida por periodicidad (eso es lo que mide el detector de diamante en subred de Zr: 0).
