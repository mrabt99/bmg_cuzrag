# Prueba de tamaño para κ: 8192 átomos

Objetivo: ver si κ depende del tamaño de celda. Méndez–Jaroszewicz usan 8192 átomos y nosotros 65536.
- **300 K:** κ de un vidrio de 8192 átomos a 1e12 K/s (protocolo del profesor, Pdamp 1.0 en el templado) contra `../cal_300K_1e12` (65536 átomos, mismo protocolo). Aísla el tamaño.
- **1000 K:** κ con 8192 átomos (`../cal_1000K_1e12_N8192`) contra `../cal_1000K_2e10` (65536). Prueba la hipótesis de propagones/cola lenta. Supone que la tasa de enfriamiento no importa en el líquido a 1000 K.

No tenemos el `in.melt` de Seba (sí se lo pedimos), así que el sistema se construye desde cero: red bcc 16×16×16 = 8192 sitios, 3891 Zr + 3891 Cu + 410 Ag con tipos asignados al azar (`set type/subset`).

## Archivos

| archivo | qué es |
|---|---|
| `in.gen8192` (v2) | construcción, fundido 1 ns a 2000 K (Pdamp 10), templado a 1e12 K/s en 17 tramos de 100 K × 100 ps (NPT iso, Tdamp 0.1, Pdamp 1.0). Escribe `restart.melt_N8192`, `restart.T1000_N8192_1e12` y `restart.quench_N8192_1e12.T300` |
| `in.relax8192` | relajación de 10 ns a 300 K y 0 Pa (NPT iso, Tdamp 0.1, Pdamp 1.0), restarts cada 1 ns. Escribe `restart.glass300K_N8192_1e12` |
| `run_gen.sbatch`, `run_relax.sbatch` | 48 tareas, `standby` (2 h y 4 h) |
| `intento1_fallido/` | log, salidas y `in.gen8192.v1` del intento que explotó |

Ambas entradas aceptan `-var fac <x>` para escalar los pasos al probar la sintaxis.

## Intento 1 (job 45120882, 2026-11-05): explotó

`in.gen8192.v1` tenía 100 ps a 2500 K y 1 ns a 2000 K, **ambos con Pdamp 1.0**, y solo un `fix halt` por T > 8000 K. Murió a los 5 min 41 s (paso 474 000, en la meseta de 2000 K) con T = 588 893 K. No llegó a escribir ningún restart, así que los jobs dependientes no corrieron (quedaron en `DependencyNeverSatisfied` y se cancelaron).

El `halt` por temperatura **solo detectó el final del colapso**. Antes hubo eventos que no vio, porque T estaba normal:

| paso | T (K) | V (Å³) | ρ (g/cm³) | Ep (eV) | qué pasa |
|---|---|---|---|---|---|
| 60 000 | 2514 | 163 156 | 6.58 | −35 889 | normal (2500 K) |
| 80 000 | 2358 | 300 283 | 3.57 | −32 523 | el volumen casi se duplica y después vuelve |
| 340 000 | 1977 | 158 616 | 6.77 | −36 248 | normal (2000 K) |
| ~350 000 | 2007 | 161 082 | 6.66 | −35 857 | Ep salta +390 eV (+48 meV/átomo) y queda ahí |
| 440 000 | 2008 | 251 920 | 4.26 | −34 520 | otra expansión |
| 449 000-473 000 | 750-7600 | ~186 000 | ~5.7 | — | T sin control |
| 474 000 | 588 893 | — | — | — | corta el `halt` |

Interpretación: es la patología conocida del potencial (`../../rate_2e10/ESTADO_PROYECTO_meltquench_1e10.md`). El EAM no tiene un núcleo repulsivo válido por debajo de ~1.2 Å. A alta T, con Pdamp chico, dos átomos pueden colapsar y quedar en un estado artificial. El salto permanente de Ep es la firma de eso. **Si la corrida no hubiera terminado explotando, habría escrito un fundido defectuoso sin que nadie lo notara.**

## Cambios en v2 (2026-11-06)

1. Sin fase a 2500 K: fundido directo a 2000 K, como el paper y el profesor.
2. Pdamp 10 en el fundido (como en 2e10/1e11), Pdamp 1.0 solo en el templado, para que el vidrio sea comparable con el de 65536 átomos a 1e12. El fundido es líquido en equilibrio, así que su Pdamp no debería afectar el vidrio final.
3. Semilla nueva (97531).
4. Detectores tempranos, cada 1000 pasos, en `in.gen8192` y en `in.relax8192`:
   - `compute pair/local dist` + `compute reduce min`, y `fix halt` si la **distancia mínima entre átomos < 1.5 Å**. En la prueba de sintaxis, a 2000 K la distancia mínima normal es ~1.97-2.0 Å, y a 300 K ~2.18 Å.
   - `fix halt` si la **densidad < 6.0 g/cm³** (normal: 6.5-7.4 según T).
   - Se mantiene el de T > 8000 K.

   La distancia mínima también se imprime en el thermo (`c_dmin`).

Riesgo que queda: los primeros tramos del templado (2000 → 1500 K) son calientes y van con Pdamp 1.0. Si vuelve a explotar, la alternativa es Pdamp 10 en todo, aceptando que en la comparación el tamaño queda mezclado con el Pdamp.

Relanzado el 2026-11-06:

| paso | jobs |
|---|---|
| generación | 45157485 |
| relajación | 45157486 |
| 300 K | NPT 45157487, NVE 45157488-90 |
| 1000 K | NPT 45157491, NVE 45157492-94 |

## Resultado (2026-11-06): κ no depende del tamaño

Todos los jobs terminaron COMPLETED. Ningún halt saltó: dmin ≥ 1.87 Å, ρ ≥ 6.73 g/cm³, sin saltos de Ep. El volumen por átomo coincide con 65536 átomos (17.565 vs 17.563 Å³ a 300 K; 18.198 vs 18.199 a 1000 K).

`compara_tamano.py` (desde `transporte/`) usa el mismo criterio que `mail_seba/kappa_vs_T.py`. Salidas: `compara_tamano.txt` y `compara_tamano.png`.

| T | N=65536 | N=8192 | diferencia |
|---|---|---|---|
| 300 K (ambos 1e12) | 0.781 ± 0.017 | 0.779 ± 0.017 | −0.002 ± 0.024 (0.1σ) |
| 1000 K (2e10 vs 1e12) | 0.864 ± 0.014 | 0.845 ± 0.016 | −0.019 ± 0.021 (0.9σ) |

Las curvas κ(t) se superponen enteras, incluido el máximo de ~2.5 W/m/K a ~70 fs a 300 K, y la cola larga es igual en los dos tamaños. El tamaño no explica la diferencia con Méndez et al.

## Recomendación general

Cualquier corrida nueva a alta T con este potencial debería llevar el `halt` por distancia mínima. El de temperatura no alcanza.
