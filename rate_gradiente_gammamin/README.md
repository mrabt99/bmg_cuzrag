# Tasa graduada hacia gamma_min (Ryltsev) — Cu-Zr-Ag BMG — Negishi/SLURM

## Qué es esto (leer antes de tocar nada)

**Corrida EXPLORATORIA.** No es parte de la comparación de tasas constantes de la
tesis (2e10/1e11/1e12 K/s, en `../rate_2e10/`, `../rate_1e11/`, `../rate_1e12_reference/`).
No tiene que salir bien: si al terminar (o al cortar antes de tiempo) el resultado es
interpretable (indicio de cristalización, o un vidrio limpio a tasa más lenta), se usa;
si explota sin recuperación, o el resultado es ambiguo, se abandona sin que esto bloquee
nada del cronograma de la tesis.

**Objetivo:** acercarse a γmin=1.5e9 K/s (Ryltsev et al., donde reportan cristalización
parcial de Cu2Zr), pero SOLO en la ventana de temperatura donde eso tiene sentido físico,
sin pagar el costo de correr todo el rango 2000→300K a una tasa ~13x más lenta que la ya
usada en la cadena 2e10.

**No se tocó ni se modificó nada de `../rate_2e10/`** — todo lo nuevo vive acá, en su
propio directorio, arrancando de una *copia* verificada (checksum) del restart de esa
cadena.

## Perfil de tasa (3 tramos, no monótono)

1. **2000→1500 K:** no se corrió nada nuevo acá. Se bifurca directamente del restart de
   la cadena 2e10 en 1500 K (`../rate_2e10/prod_2e10/restart.quench_N65536_2e10.seg4` —
   confirmado por el propio input `in.quench_2e10_seg04.lmp`, que dice "1600 K -> 1500 K"
   y escribe ese restart; el siguiente, `seg05.lmp`, arranca desde ahí hacia 1400K). Esta
   zona ya está pagada y es irrelevante para la estructura final. Copia local verificada
   bit-a-bit (sha256 idéntico): `restart.start_1500K_from_2e10_seg4`.

2. **1500→700 K (tramos 1-8, `in.quench_grad_seg01..08.lmp`):** ventana que importa —
   donde Ryltsev ve el efecto de la tasa y donde puede aparecer cristalización. Tasa
   interpolada logarítmicamente entre 2e10 K/s (tramo 1) y 1.5e9 K/s = γmin (tramo 8),
   factor **1.4478** entre tramos consecutivos (~1.45 esperado). Mismo potencial EAM sin
   parchear, mismo ensemble NPT Pdamp=10.0/Tdamp=0.1, dt=1fs — todo derivado de
   `../rate_2e10/prod_2e10/in.quench_2e10_seg04.lmp`, no reescrito desde cero.

3. **700→300 K (tramos 9-12, `in.quench_grad_seg09..12.lmp`):** tasa rápida, **1e11 K/s**
   (elegida por el usuario sobre la alternativa 2e10, para no pagar walltime extra en una
   zona no-ergódica donde la tasa ya no afecta la estructura final). Derivados de
   `../rate_1e11/prod_1e11/in.quench_1e11_seg09.lmp`.

## Tabla real de tasas y duraciones

| tramo | T inicio→fin | tasa (K/s) | pasos | duración | walltime pedido (sbatch) |
|---|---|---|---|---|---|
| 1 | 1500→1400 | 2.000e10 | 5,000,000 | 5.00 ns | 10:42:00 |
| 2 | 1400→1300 | 1.381e10 | 7,238,949 | 7.24 ns | 15:15:00 |
| 3 | 1300→1200 | 9.542e9 | 10,480,478 | 10.48 ns | 21:51:00 |
| 4 | 1200→1100 | 6.590e9 | 15,173,529 | 15.17 ns | 1-07:25:00 |
| 5 | 1100→1000 | 4.552e9 | 21,968,082 | 21.97 ns | 1-21:15:00 |
| 6 | 1000→900 | 3.144e9 | 31,805,166 | 31.81 ns | 2-17:18:00 |
| 7 | 900→800 | 2.172e9 | 46,047,198 | 46.05 ns | 3-22:18:00 |
| 8 | 800→700 (=γmin) | 1.500e9 | 66,666,667 | 66.67 ns | 5-16:19:00 |
| 9 | 700→600 | 1.000e11 | 1,000,000 | 1.00 ns | 02:33:00 |
| 10 | 600→500 | 1.000e11 | 1,000,000 | 1.00 ns | 02:33:00 |
| 11 | 500→400 | 1.000e11 | 1,000,000 | 1.00 ns | 02:33:00 |
| 12 | 400→300 | 1.000e11 | 1,000,000 | 1.00 ns | 02:33:00 |

Total ventana graduada: ~204.4M pasos (~204 ns). Walltime pedido = ~2.2x el tiempo
estimado a ~300 pasos/s (48 ranks), mismo criterio de margen que `rate_2e10`/`rate_1e11`.
QOS `cpu-qpart` en Negishi permite hasta 14 días por job, así que ningún tramo (ni el 8,
el más largo) necesitó partirse en sub-jobs.

## Cinturones de seguridad implementados

- `fix halt_T` (T > 8000K) y `fix halt_P` (|P| > 1e5 bar), igual que en `rate_2e10`/`rate_1e11`.
- **Recuperación ante explosión — MODIFICADA respecto a la cadena 2e10:** en la zona
  1500→700K NO se usa hold isotérmico (un hold cerca de Tg es envejecimiento extra en la
  variable que se está midiendo). Si un tramo explota: volver al restart de ENTRADA de
  ESE tramo y relanzar con `velocity all create <T> <nueva_seed> mom yes rot yes dist
  gaussian` justo antes del `fix qnpt`, sin espera adicional a esa T. Documentado en el
  encabezado de cada `in.quench_grad_segNN.lmp` (tramos 1-8). Si el mismo tramo explota
  más de 2-3 veces: no insistir, es señal válida (puede ser la cercanía a la zona de
  cristalización/inestabilidad de γmin) — reportar y cortar ahí.
- `check_crystallization.py`: monitor de caída conjunta de energía potencial + salto de
  densidad por encima de la tendencia suave de enfriamiento (umbral configurable, default
  8σ robusto — mismo criterio usado ad hoc en el análisis de `rate_1e11`). Uso:
  `python3 check_crystallization.py quench_N65536_gradmin.log.dat.seg0*`. Si señala algo,
  no seguir la cadena más allá de ese tramo sin avisar.
- `chain_submit.sh` soporta `--from N` y `--upto N` para poder pausar la cadena entre
  tramos si hiciera falta revisar antes de seguir a los tramos más lentos.
- Test de sanity (3000 pasos, sin NaN/Inf/lost atoms) corrido ANTES de lanzar la cadena
  real sobre las 4 configuraciones representativas (tramo 1 = tasa rápida/Nfreq chico,
  tramo 8 = tasa más lenta/Nfreq grande, tramos 9 y 12 = cola 1e11) — ver `test/`.
  Resultado: limpio en las 4 (job 43306374).

## Decisión del usuario sobre el modo de ejecución (2026-09-18)

Se lanzó la cadena completa de una sola vez con `--dependency=afterok` (los 12 tramos,
job IDs 43306643-43306655), confiando en `fix halt` como cinturón de seguridad, en vez de
un esquema de lanzamiento tramo-por-tramo con revisión automática de cristalización entre
cada uno. **La revisión de `check_crystallization.py` sobre cada tramo recién terminado
queda a cargo del usuario**, no es automática en esta cadena.

## Aviso operativo: mantenimiento de Negishi

Hay una reserva de mantenimiento de todo el cluster (`Sept_Power_Maint`) que arranca
2026-09-19T00:00:01. El tramo 1 corre y termina antes de esa hora sin problema, pero los
tramos 2-7 (cuyo walltime pedido los haría cruzar esa ventana si arrancaran hoy) van a
quedar en cola esperando a que Negishi vuelva — esto es normal, la cadena `afterok` no se
rompe por esto, simplemente se pausa. El `EndTime` que reporta `scontrol show reservation`
(2027-09-19) es casi seguro un placeholder administrativo de RCAC, no la duración real del
corte.

## Estructura de archivos

```
ZrCuAg.lammps.eam                        <- potencial, copiado de rate_2e10 (sin modificar)
restart.start_1500K_from_2e10_seg4       <- copia verificada (sha256) del restart 1500K
                                              de rate_2e10/prod_2e10/restart.quench_N65536_2e10.seg4
in.quench_grad_seg01.lmp .. seg08.lmp    <- ventana graduada 1500->700K
in.quench_grad_seg09.lmp .. seg12.lmp    <- cola rapida 700->300K a 1e11 K/s
slurm_seg01.sbatch .. seg12.sbatch       <- un sbatch por tramo
chain_submit.sh                          <- encadena con afterok, soporta --from/--upto
check_crystallization.py                 <- monitor PE+densidad para la ventana graduada
test/                                    <- sanity test de 3000 pasos (job 43306374), limpio
```

## Al terminar (o al cortar antes de tiempo)

Pendiente completar en este README: tabla real de qué tramos corrieron limpios a la
primera vs. cuántos necesitaron reintento con nueva semilla, en qué tramo (si alguno)
apareció señal de `check_crystallization.py`, y gráfico de volumen/energía por átomo vs.
T de punta a punta (mismo pipeline que `../structural_analysis/`).
