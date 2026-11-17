# rate_1e12: vidrio a 1e12 K/s con el MISMO protocolo que 2e10 y 1e11 (16/11/2026)

Esta carpeta reemplaza a `rate_1e12_reference/`, la referencia del profesor, que **no se borra**. Objetivo: que las tres tasas tengan la misma metodología y que la tasa sea la única variable.

## Qué cambia respecto de `rate_1e12_reference/`

| | `rate_1e12_reference/` (profesor) | `rate_1e12/` (esta) = 2e10 = 1e11 |
|---|---|---|
| punto de partida | `restart.melt_N65536` (1 ns NPT a 2000 K, Pdamp 1) | `restart.equil2000K_N65536_1e10` (+2 ns a 2000 K, semilla "seed2") |
| 2000 → 1900 K | seg1 de su cadena | barrido `sweep_1e12/` (ya hecho el 05/09, job 42877510, Pdamp 10) |
| Pdamp en el temple | 1 ps | **10 ps** |
| tramos | 17 de 100 K (LAMMPS 2018, log crudo) | 1 barrido + 16 tramos de 100 K (fix ave/time, mismo formato que 1e11) |
| relajación a 300 K | 10 ns (Pdamp 1) | 10 ns, NPT Pdamp 10 (`in.relax300K_1e12` = `in.relax300K_1e11` con la tasa cambiada) |

## Contenido

- `ZrCuAg.lammps.eam`: copia de `rate_1e11/`, con el mismo sha256 que `rate_2e10/`.
- `sweep_1e12/`: copia de `rate_2e10/sweep_1e12/`. El restart tiene el mismo sha256 que el original (fed403b9…).
- `prod_1e12/`:
  - `in.quench_1e12_seg01..16.lmp`, generados desde `rate_1e11/prod_1e11/in.quench_1e11_seg02.lmp` cambiando solo la tasa: 1e5 pasos por tramo (0.1 ns, 100 K), dump cada 1e4 pasos (10 cuadros por tramo, cada 10 K). El tramo 1 lee `../sweep_1e12/restart.sweep_1e12_1900K`.
  - `in.relax300K_1e12`.
  - `slurm_ramp.sbatch`: los 16 tramos en un solo job, `standby`, 4 h. Si un tramo falla, se retoma con `sbatch --export=ALL,FROM=N slurm_ramp.sbatch`.
  - `slurm_relax.sbatch`: QOS normal, 14 h. La de 1e11 tardó 9 h.
- Salidas: las mismas que 1e11 (`restart.quench_N65536_1e12.segN`, `quench_N65536_1e12.log.dat.segN`, `.dump.segN`, `restart.glass300K_N65536_1e12`, `relax300K_N65536_1e12.log.dat`).

## Lanzamiento (16/11/2026, a pedido del usuario)

- Rampa: job 46018926.
- Relajación: job 46018927 (`afterok`).
- Costo estimado: ~1.7 h la rampa y ~9 h la relajación.
- **El primer intento falló en el tramo 3** (16/11/2026). Fue un error al generar los inputs: `in.quench_1e12_seg03.lmp` leía `seg3` en vez de `seg2`. Los tramos 1 y 2 terminaron bien (el 2 cerró a T = 1701.6 K, ρ = 6.895 g/cc). Se corrigió el input, se revisó la cadena de los 16 y se canceló el relax colgado (46018927, `DependencyNeverSatisfied`). Relanzado desde el tramo 3 (`FROM=3`): rampa **46020228**, relax **46020229** (`afterok`). Los logs del primer job llevan el sufijo 46018926.

## Resultado (17/11/2026)

Rampa (job 46020228, 1 h 24 min) y relajación (job 46020229, 8 h 40 min) COMPLETED, sin warnings. Los 16 tramos tienen 100 líneas de log cada uno y están todos los restarts. Vidrio final: T = 300.0 K, Ep = −4.83018 eV/at, ρ = 7.4615 g/cc, Lx = 104.796 Å.

Thermo en `structural_analysis/1e12u/thermo_1e12u.csv` (extract_thermo.py, avetime: sweep + 16 tramos + relax, 3600 filas). Tg en `1e12u/tg_1e12u.py` → `Tg_1e12u.csv`, con compute_tg.py sin cambios:

| | Tg(V) | Tg(Ep) | media | Ep final (eV/at) | ρ final |
|---|---|---|---|---|---|
| 1e12 unificada | 698.0 | 698.0 | 698.0 | −4.83018 | 7.4615 |
| 1e12 referencia vieja | 703.2 | 697.9 | 700.5 | −4.82908 (a T = 302.3 K) | 7.4593 |

El cambio de protocolo (Pdamp 1 → 10 y punto de partida) mueve muy poco las cosas. El vidrio nuevo queda apenas más relajado: Ep ≈ −0.8 meV/at corrigiendo los 2.3 K de diferencia, ρ +0.03 %, Tg −2.5 K. El orden entre tasas no cambia.

Manifest de 17 restarts: `structural_analysis/1e12u/manifest.txt` (copia en `etapa4/manifests/1e12u.txt`). La etiqueta de los análisis nuevos es **1e12u**.

## Hecho el 17/11/2026: estructura (etapa A)

Ver `structural_analysis/1e12u/README.md`: dumps, BOO, Voronoi, RDF/BADF y tablas. Las conclusiones no cambian (Cu-ico BOO 2e10/1e12: 1.61 contra 1.74 con la vieja). La bitácora (secciones 2-6) y `estilo/bmg_estilo.py` (TG_PROPIAS 1e12 = 698.0) ya usan esta corrida.

## Pendiente (no hecho)

Los puntos 1-3 de abajo ya están hechos; quedan el 4 (mecánica, Cij y transporte: etapas B y C, que necesitan el OK del usuario) y el 5 en lo que respecta a las secciones 7-9 de la bitácora.


Nada de lo de abajo se hizo todavía. Todos los análisis existentes de "1e12" usan la referencia vieja.

1. Chequear la integridad de los logs y que T, Ep, V y ρ a 300 K sean razonables. Comparar con la referencia vieja: ¿cuánto cambia por Pdamp 1 → 10 y por el punto de partida?
2. Thermo con `extract_thermo.py --style avetime` y Tg con `compute_tg.py`, en carpetas nuevas.
3. Manifest de 17 restarts (`sweep_1e12/restart.sweep_1e12_1900K` + `prod_1e12/seg1..15` + `glass300K`), y con él dumps, BOO, Voronoi y RDF/BADF.
4. Mecánica, Cij y transporte con el vidrio nuevo.
5. Bitácora: reemplazar 1e12 en todas las secciones.
