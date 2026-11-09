# κ por NEMD (método directo con `fix ehex`) contra Green-Kubo, vidrio a 300 K

**Objetivo:** medir κ con un método independiente de Green-Kubo. Si NEMD da ≈ 0.78 W/m/K, confirma la meseta de GK y descarta el valor de Méndez et al. (2.43). Ese valor coincide con el máximo de κ(t), a ~70 fs.

## Sistema y referencia

- **Vidrio:** el de 8192 átomos templado a 1e12 y relajado 10 ns a 300 K (`../size_8192/restart.glass300K_N8192_1e12`). Se lleva al volumen medio de `../cal_300K_1e12_N8192` (143 891.9 Å³, lado 52.40 Å) y se replica n veces en z.
- **Referencia GK, mismo vidrio:** κ = 0.779 ± 0.017 W/m/K. Con 65536 átomos al mismo protocolo da 0.781 ± 0.017, así que GK no muestra efecto de tamaño (`../size_8192/compara_tamano.txt`).

| n | L_z (Å) | átomos | tramos de 1 ns | transitorio descartado |
|---|---|---|---|---|
| 4 | 210 | 32 768 | 2 | 0.5 ns |
| 8 | 419 | 65 536 | 3 | 1 ns |
| 16 | 838 | 131 072 | 5 | 3 ns |

## Método

- **Equilibrado:** 50 ps de NVT a 300 K. Después, NVE con `fix ehex`: +Q = 0.5 eV/ps en la región caliente y −Q en la fría.
- **Reservorios:** dos losas de 10 Å, separadas L/2. La fría está en [z_lo, z_lo + 10 Å] y la caliente en [z_lo + L/2, z_lo + L/2 + 10 Å]. Elegí `ehex` y no Müller-Plathe porque la aleación tiene tres masas distintas; con `ehex` el flujo queda impuesto exacto y la deriva de energía es mínima.
- **Flujo:** se reparte en dos mitades, así que J = Q/(2A), con A = 52.40² Å². En unidades prácticas, κ[W/m/K] = 801.09 · Q[eV/ps] / (A[Å²] · |dT/dz|[K/Å]).
- **Gradiente esperado** (κ ≈ 0.78): 0.19 K/Å. La diferencia de T entre reservorios es ~20 K (n=4), ~40 K (n=8) y ~80 K (n=16; T entre ~260 y 340 K, todo vidrio y lejos de Tg).
- **Tiempo de relajación:** la difusividad térmica es α ≈ κ/(3 n k_B) ≈ 33 Å²/ps, y el tiempo del modo lento, (L/2)²/(π²α), da ≈ 35, 135 y 540 ps.
- **Perfil:** T(z) en bins de 5 Å, promediado en bloques de 100 ps (`fix ave/chunk ... temp`).
- **Ajuste:** recta en el 60 % central de cada mitad. El sistemático sale de variar ese rango a 50 % y 70 %.
- **Error estadístico:** dispersión entre superbloques de 500 ps del régimen estacionario.
- **Tamaño:** 1/κ contra 1/L para extrapolar κ_∞. En un amorfo la dependencia con L debería ser débil; una subida con L indicaría contribución de propagones.
- **Halts:** distancia mínima < 1.5 Å, densidad < 6.0 g/cm³ y T > 2000 K.

## Archivos

| archivo | qué es |
|---|---|
| `in.nemd` | entrada por tramos de 1 ns, encadenada con restarts (`-var n`, `-var part`, `-var seed`, `-var Q`, `-var fac` para pruebas) |
| `run_nemd.sbatch` | un tramo; variables de entorno N, PART, SEED (`standby`, 48 tareas, 4 h) |
| `lanzar_nemd.sh` | encola las tres cadenas con `afterok` (lanzar solo con go-ahead) |
| `analiza_nemd.py` | perfil promedio, κ por longitud, κ por bloque contra t, extrapolación 1/L; escribe `nemd.txt` y `nemd.png` |
| `prueba_sintaxis/` | prueba corta en el login (n=2, `fac 0.01`, partes 1 y 2) |

## Prueba de sintaxis (2026-11-06)

n=2 (16 384 átomos, Lz = 104.8 Å), `fac 0.01`, 8 procesos en el login:
- Las partes 1 y 2 terminan con rc=0. La parte 2 lee el restart y sigue el reloj (pasos 10 000 → 20 000).
- E_total se conserva exacta (−78 504.4 eV en las dos partes).
- A los 10 ps el perfil ya muestra el gradiente: ~292 K junto a la losa fría y ~306-311 K junto a la caliente.
- `analiza_nemd.py` lee bien los bloques y la caja. Con 20 ps el κ no significa nada.

**Costo estimado** con 48 tareas, a partir del GK de 8192 átomos (3.4 × 10⁻⁶ CPU·s por átomo·paso; en el login fue ~2.7 veces más lento):
- n=4: ~0.65 h/ns
- n=8: ~1.3 h/ns
- n=16: ~2.6 h/ns

Total ~18 nodo·h. La cadena de n=16 son 5 tramos de ~2.6 h en serie (~13 h).

Salidas por tramo: `L<n>_p<k>_Tprof.dat`, `L<n>_p<k>_thermo.dat`, `L<n>.part<k>.restart`, `log.L<n>_p<k>.<job>`.

## Lanzado (2026-11-06)

| n | jobs (tramos encadenados con afterok) |
|---|---|
| 4 | 45171535-36 |
| 8 | 45171537-39 |
| 16 | 45171540-44 |

## Resultado (2026-11-07)

`analiza_nemd.py 4 8 16` → `nemd.txt`, `nemd.png`. Todos los jobs COMPLETED; ningún halt saltó.

| L (Å) | κ (W/m/K) | ΔT entre reservorios |
|---|---|---|
| 210 | 0.695 ± 0.005 ± 0.010 | 21 K |
| 419 | 0.709 ± 0.007 ± 0.005 | 42 K |
| 838 | 0.722 ± 0.003 ± 0.005 | 85 K |
| ∞ (1/κ vs 1/L) | **0.731 ± 0.009** | — |

- Los errores son estadístico (superbloques) y sistemático (rango de ajuste); χ² de la extrapolación = 0.14.
- GK del mismo vidrio: 0.779 ± 0.017 (± 0.019 de ventana). La diferencia es 0.048 (6 %): 1.8σ si se incluye el sistemático de GK.
- **NEMD confirma GK y descarta el 2.43 de Méndez et al.**

## Prueba de linealidad y revisión de errores (2026-11-09)

- `L8q025` (Q = 0.25, 5 × 1 ns): κ = 0.708 ± 0.032, igual a Q = 0.5. Ver `compara_Q.txt/png`.
- Los errores por superbloques de la tabla de arriba están probablemente subestimados (bloques correlacionados, modos térmicos lentos).
  Con el ruido estimado: ±0.027/0.017/0.009 y **κ∞ = 0.732 ± 0.017**; la dependencia con L no es significativa.
- `lanzar_nemd_L8ext.sh`: extensión de L8 (tramos 4-6, Q = 0.5) para medir ese ruido.
