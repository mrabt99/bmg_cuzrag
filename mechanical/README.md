# mechanical — ensayos de tracción sobre los vidrios Zr47.5 Cu47.5 Ag5 (N = 65536)

Esta carpeta es nueva (2026-10-07). Acá vive todo lo mecánico: el ensayo **piloto** que valida el procedimiento y, más adelante,
la matriz completa (varias tasas de enfriamiento, temperaturas y tasas de deformación). **No toca** `rate_*`, `structural_analysis/` ni `qgrad*`:
de ahí solo se LEEN los restarts de los vidrios y el potencial.

## Contenido

```
mechanical/
├── README.md                      este archivo
└── pilot_tension/                 ensayo piloto: 2e10, 300 K, 1e9 s^-1, tracción en x, hasta ε = 0.2
    ├── in.tension_pilot           input de LAMMPS (reusable para la matriz con -var)
    ├── slurm_tension_pilot.sbatch script de envío al cluster
    ├── smoke/                     prueba corta (6000 pasos) que validó el input; no es parte de los resultados
    │   └── intento1_fallido/      primer intento (falló a los 11 s, ver "Historial")
    └── (salidas del piloto: log.*, *.curve.dat, *.dump, *.restart, ver abajo)
```

## Qué hace el ensayo (resumen)

1. **Etapa 1, mantenimiento** (50 ps = 50 000 pasos). `fix npt` isotrópico a T y presión 0, Pdamp = 1.0 ps. Sirve para que el vidrio se asiente
   en este input antes de estirarlo. Al terminar se escribe `<tag>.ref.restart` (estado de referencia) y se **congelan** `L0 = lx`, `Ly0 = ly`, `Lz0 = lz`.
2. **Etapa 2, tracción** (200 000 pasos = ε de 0 a 0.2). `fix deform` estira x a velocidad constante (`erate`, `units box`, `remap x`);
   `fix npt` solo en y y z con presión 0 (la caja se contrae libre: efecto Poisson). Eso es tracción uniaxial.

Se usan las velocidades del restart (no se recrean). `timestep 0.001` (1 fs) **no se cambia**: hace falta para el análisis de fonones posterior.
El potencial y `pair_coeff * * ... Zr Cu Ag` son los mismos del temple, para que el mapeo de tipos (1 = Zr, 2 = Cu, 3 = Ag) coincida con el restart.

## Parámetros del input (`-var nombre valor`)

| variable | por defecto | qué es |
|---|---|---|
| `sample` | 2e10 | tasa de enfriamiento del vidrio: elige el restart (2e10, 1e11 o 1e12) |
| `T` | 300 | temperatura (K) |
| `erate` | 1e9 | tasa de deformación en **1/s** (se convierte a 1/ps: 1e9 s^-1 = 1e-3 ps^-1) |
| `pdamp` | 1.0 | Pdamp del barostato (ps) |
| `emax` | 0.2 | deformación final |
| `nhold` | 50000 | pasos de mantenimiento (etapa 1) |

Todo lo demás se deriva de `erate` (deformación por paso = erate·dt = 1e-6 con 1e9 s^-1):
pasos de la etapa 2 = emax/(erate·dt) = 200 000; un punto de curva cada Δε = 0.001 (1000 pasos), un dump cada Δε = 0.005 (5000 pasos),
un restart cada Δε = 0.05 (50 000 pasos). Válido para `erate ≤ 1e10`. Dirección fija en x por ahora.

Los restarts de origen (solo lectura):
- 2e10: `rate_2e10/prod_2e10/restart.glass300K_N65536_2e10`
- 1e11: `rate_1e11/prod_1e11/restart.glass300K_N65536_1e11`
- 1e12: `rate_1e12_reference/restart.glass300K_N65536_1e12`
- potencial: `rate_2e10/ZrCuAg.lammps.eam` (idéntico, mismo md5, al de qgrad_ext)

## Salidas (con `tag = tens_<sample>_<T>K_<erate>_x`, p. ej. `tens_2e10_300K_1e9_x`)

| archivo | contenido |
|---|---|
| `<tag>.curve.dat` | curva σ-ε: **una línea por Δε = 0.001**, promediada dentro del intervalo (100 muestras). Columnas: `step strain sxx syy szz T PE(eV/átomo) lx ly lz vol`. Tensiones en **GPa** (σ = −p/10000, tracción positiva); `strain` = (lx − L0)/L0 (ingeniería), también promediado. La primera línea tiene ε ≈ 0.0005 y no mezcla datos de la etapa 1. |
| `<tag>.dump` | `id type x y z vx vy vz`, ordenado por id, cada Δε = 0.005, **incluye ε = 0** (41 cuadros). Los IDs son necesarios para la deformación atómica en OVITO. |
| `<tag>.ref.restart` | estado de referencia (fin del mantenimiento, ε = 0) |
| `<tag>.<paso>.restart` | restarts cada Δε = 0.05: pasos 100000, 150000, 200000, 250000 (ε = 0.05, 0.10, 0.15, 0.20) con nhold = 50000 |
| `log.<tag>.<jobid>` | log de LAMMPS (imprime `L0`, `Ly0`, `Lz0` y las líneas `###`; thermo cada 1000 pasos en ambas etapas) |
| `slurm_tens_pilot.<jobid>.out/.err` | salida estándar y de error de Slurm (el `.err` solo trae el aviso habitual de módulos/xalt) |

Nota: el esfuerzo exacto a ε = 0 no está en la curva (su primer punto promedia 0–0.001); está en el thermo del paso `nhold`.

## Medidas de seguridad (`fix halt`, cada 1000 pasos)

- `temp > 1500 K` → aborta.
- `ly > 1.1·Ly0` o `lz > 1.1·Lz0` → aborta (en tracción ly y lz solo pueden achicarse; que crezcan indica un problema).
- No hay halt por densidad: en tracción la densidad baja legítimamente.

## Cómo correrlo

```bash
cd /scratch/negishi/aaabt/bmg_cuzrag/mechanical/pilot_tension
sbatch slurm_tension_pilot.sbatch        # 48 ranks, standby, 1 h
squeue -u $USER                           # estado
```
Medido: ≈ 330 pasos/s con 48 ranks → 250 000 pasos ≈ 13–20 min. `sbatch`: `-A strachan -p cpu -q standby`, `srun --mpi=pmi2`
(pmix dejó de andar tras el mantenimiento de septiembre), módulos `gcc/12.2.0 openmpi/4.1.4 lammps/20220623`.
Regla del grupo: jobs de menos de 4 h van con `-q standby`. **No lanzar nada sin confirmar** (regla de trabajo).

Para otra muestra/temperatura/tasa, editar las tres variables `SAMPLE`, `TEMP`, `ERATE` del sbatch (o pasar `-var` directo a `lmp`) y correr desde una carpeta nueva
para no mezclar salidas (las salidas se escriben en el directorio desde donde se hace `sbatch`).

## Historial y verificación

- **Prueba corta** (`smoke/`, `-var nhold 1000 -var emax 0.005`, 6000 pasos, 23 s): sin warnings; ε llegó a 0.005; T entre 298 y 301 K;
  syy y szz ≈ 0; Ly y Lz se achican (Poisson ≈ 0.38); dump con 2 cuadros (ε = 0 y 0.005); curva con 5 líneas, la primera en ε ≈ 0.0005.
- **Primer intento fallido** (`smoke/intento1_fallido/`, job 43896653): `dump` no admite una variable como frecuencia en este LAMMPS (20220623).
  Se corrigió usando `dump_modify ... every v_dnext` (y `restart v_rnext`), que alinean dumps y restarts con el inicio de la etapa 2 aunque `nhold` no sea múltiplo de 5000/50000.
- **Piloto** (job 43896715): lanzado 2026-10-07 tras la prueba corta. El primer restart de la tracción (paso 100000) se escribió bien, así que `restart v_rnext` funciona.
  Cuando termine, ver `squeue`/`sacct -j 43896715` y el final del log (`### TRACCION TERMINADA`).

## Resultados del piloto (job 43896715, COMPLETED, 13 min 56 s, sin warnings; 41 cuadros de dump, 200 puntos de curva, ε final = 0.2000)

Análisis con `analyze_curve.py` (salidas en `pilot_tension/analysis/`: `curva_sigma_eps.png`, `curva_resumen.txt`):

| magnitud | valor |
|---|---|
| E (ajuste lineal, ε ≤ 0.01, 10 puntos, ordenada libre) | 54.6 ± 0.5 GPa (por el origen: 60.4 GPa) |
| σmax | 2.27 GPa en ε = 0.0615 (con suavizado de 5 puntos: 2.26 GPa en ε = 0.0635) |
| σ de flujo (media 0.15 ≤ ε ≤ 0.20) | 1.30 ± 0.04 GPa |
| Δσ = σmax − σflujo | 0.97 GPa |
| verificación | \|σyy\|, \|σzz\| ≤ 0.005 GPa; T = 300.04 ± 0.26 K; volumen +0.72 % al final; ly y lz bajan a 95.7 y 96.3 Å (de 104.75) |

Los "±" son el error estándar del ajuste y el desvío dentro de la ventana: ruido térmico de UNA realización, no variación entre muestras.
Coeficiente de Poisson inicial ≈ 0.41 (de la contracción lateral con ε ≤ 0.01) → G ≈ 19 GPa y K ≈ 97 GPa si el material fuera isótropo elástico.

Deformación atómica de von Mises con `von_mises_strain.py` (OVITO, cuadro de ε = 0.20 respecto de ε = 0; `analysis/eta_*`):
- Cutoff = primer mínimo de la g(r) total del vidrio 2e10 a 300 K = **3.75 Å** (primer pico en 2.81 Å; se calcula de `structural_analysis/etapa4/rdf/rdf_all.csv` suavizando 5 bins, porque viene de un solo snapshot; sin suavizar da 3.83 Å).
- OVITO: `AtomicStrainModifier`, referencia = cuadro 0, `affine_mapping = ToReference`, `eliminate_cell_deformation = False` (η es la deformación TOTAL; con `Off` da lo mismo, y con `eliminate_cell_deformation = True` la media baja a 0.23).
- η: media 0.290, mediana 0.266, desvío 0.149, p95 0.57, p99 0.75, máximo 1.61. Verificación: la media coincide con el valor homogéneo esperado de Green-Lagrange con las dimensiones de la caja (0.300); el cuadro 0 contra sí mismo da η = 0.
- Referencia de ruido: de ε = 0 a 0.005 la mediana de η es 0.04 y p99 = 0.09 (curva gris del histograma).
- 9.0 % de los átomos tiene η > 0.5 y 0.09 % tiene η > 1.0.
- En los cortes finos (10 Å, vistas x-y y x-z, escala de color saturada en 1.0) la deformación es heterogénea, con parches de η alto; **a simple vista no se ve una banda de corte localizada** a ε = 0.20 (ver lectura).

Para repetir:
```bash
cd pilot_tension
/scratch/negishi/aaabt/venv_structural/bin/python analyze_curve.py            # numpy + matplotlib
module load anaconda/2024.02-py311; python von_mises_strain.py               # necesita OVITO (el venv no lo tiene)
```
`eta_eps0.20.npz` guarda `id` y `eta` por átomo (útil para correlacionar con la estructura previa).

## Pendiente

1. Matriz completa: reusar `in.tension_pilot` con `-var sample/T/erate` (y carpeta nueva por caso).
2. Decidir si el piloto basta como validación o conviene repetirlo con otra semilla de velocidades / otro vidrio, para tener una barra de error entre muestras.

## Cosas a revisar en los resultados del piloto

- E ≈ 55 GPa con ν ≈ 0.41 (G ≈ 19 GPa) me parece **bajo/blando** para un vidrio Zr-Cu (de memoria, sin verificar contra bibliografía: E ≈ 70–90 GPa, ν ≈ 0.36–0.38, G ≈ 30 GPa).
  El ajuste sobre ε ≤ 0.01 ya incluye algo de no linealidad, el ensayo es isotérmico (E isotérmico < adiabático) y la composición y el potencial no son los de esas referencias,
  pero conviene contrastarlo (p. ej. E con ε ≤ 0.005, o constantes elásticas por deformaciones pequeñas) antes de usar los valores absolutos de la matriz.
- σmax = 2.27 GPa con un pico ancho en ε ≈ 0.06 y reblandecimiento gradual hasta 1.3 GPa, sin caída brusca ni banda de corte clara a ε = 0.2: puede ser lo esperable a 1e9 s⁻¹ y 300 K en una caja de 10 nm,
  pero el criterio con el que se juzgó (inspección de dos cortes) es cualitativo.
- Pdamp = 1.0 ps es mucho más corto que el 10 ps del temple (que se había subido porque el barostato explotaba durante el enfriamiento del líquido);
  a 300 K en un vidrio no debería dar problemas, y el halt de ly/lz y de T lo detectaría.
- `remap x` mueve posiciones pero no velocidades: a 1e9 s^-1 el flujo de caja es ≈ 10 m/s frente a cientos de m/s de agitación térmica, efecto pequeño.
- Una sola muestra (un vidrio por tasa) y una sola realización de la tracción: no hay barras de error entre muestras.

## Ampliaciones (2026-10-07)

- **`scripts/localization.py`**: localización cuantitativa (η de von Mises INCREMENTAL en 10 ventanas consecutivas de Δε = 0.02 entre 0 y 0.20, referencia en el inicio de cada ventana;
  ψ = desvío estándar de η entre átomos; figura σ-ε arriba y ψ abajo con eje ε compartido; render de la caja entera con η > umbral (0.2, `--thr`) en ε 0.06→0.10 y 0.16→0.20, vistas a lo largo de y y z;
  `.npz` con id y η por ventana). Uso: `python scripts/localization.py --dir pilot_tension --rate 2e10` (con el python de anaconda, necesita OVITO). Salidas en `<dir>/analysis/localization/`.
- **`pilot_tension_1e12/`**: mismo piloto con el vidrio 1e12 (job 43898161, COMPLETED, 15 min; input idéntico por md5 al de `pilot_tension/`). El log trae un solo warning, `Old restart file format revision. Switching to compatibility mode` (el restart de 1e12 es de un formato más viejo; no afecta).
- **Comparación 2e10 vs 1e12** (una realización cada una): E = 54.6 vs 50.6 GPa; σmax = 2.27 (ε 0.0615) vs 1.81 GPa (ε 0.0585); σflujo = 1.30 vs 1.22 GPa; Δσ = 0.97 vs 0.59 GPa. ψ: ≈ 0.019–0.021 en la primera ventana (casi elástica) y se estabiliza en ≈ 0.050 desde ε ≈ 0.08 en ambas;
  el 1e12 sube antes (ψ = 0.028 vs 0.021 en 0.02–0.04; 0.037 vs 0.032 en 0.04–0.06). No hay banda localizada visible en los renders de ninguno (≈ 12–13 % de átomos con η > 0.2 en ventanas de 0.04).
- **`elastic_2e10/`** (inputs escritos, SIN correr): constantes elásticas del vidrio 2e10 adaptando `examples/ELASTIC_T/DEFORMATION` (300 K: `300K/`) y `examples/ELASTIC` (0 K: `0K/`),
  ver `elastic_2e10/*/init.mod`. Los ejemplos originales de LAMMPS (stable_23Jun2022) están en `pilot_tension/smoke/lammps_src/examples/`.
  **Pruebas cortas hechas (2026-10-07)** en `300K/smoke*/` y `0K/smoke/` (solo columnas 1 y 4; sin significado estadístico, validan el flujo): a 300 K C11 ≈ 131, C12 ≈ 92, C44 ≈ 19.4–19.8 GPa (→ K ≈ 105, G ≈ 19.5, E ≈ 55, ν ≈ 0.41: coincide con el piloto de tracción);
  a 0 K (δ = 1e-4) C11 = 142.8, C12 = 96.6, C44 = 23.3 GPa (→ E ≈ 65, ν ≈ 0.40); presiones de ± antisimétricas a 0.1 bar. La prueba corta reveló un error en el promedio (ventanas de `fix ave/time` desalineadas → ceros en una corrida real):
  corregido con ventanas de 1 ps y `ave running`, validado con conteos desalineados a propósito (`300K/smoke2/`). Medidas de costo: NVT con caja triclínica ≈ 245 pasos/s (→ 300 K ≈ 110 min, walltime 3 h); preparación 0 K ≈ 30 min (hecha una sola vez en `0K/prep/`), cada columna a 0 K tarda segundos.
  Post-proceso: `elastic_2e10/elastic_post.py` (matriz Cij, simetría, isotropía, K/G/E/ν de Voigt-Reuss-Hill, v_L, v_T, ΘD de Anderson, errores por bloques a 300 K); probado con una matriz isótropa sintética (K = 100, G = 20 → E = 56.25, ν = 0.4063, ΘD = 212 K, correcto). Corridas completas: pendientes de OK.

## Constantes elásticas del vidrio 2e10 — resultados (jobs 43899559 y 43899560, ambos COMPLETED, sin warnings; `elastic_2e10/`, post-proceso `elastic_post.py`)

| | 300 K (NVT, δ = ±0.005, 100 ps por estado) | 0 K (δ = ±1e-4) | 0 K (δ = ±5e-4) |
|---|---|---|---|
| C11 / C12 / C44 (promedios de C11,C22,C33 / C12,C13,C23 / C44,C55,C66) [GPa] | 131.2 / 92.4 / 19.4 | 142.7 / 96.6 / 23.1 | 142.5 / 96.7 / 23.0 |
| K, G (Hill) [GPa] | 105.3, 19.39 | 112.0, 23.09 | 112.0, 23.00 |
| E [GPa], ν | 54.8, 0.413 | 64.8, 0.4035 | 64.6, 0.404 |
| v_L, v_T [m/s] | 4190, 1611 | 4347, 1748 | 4345, 1745 |
| ΘD de Anderson [K] | 209.3 | 227.7 | 227.3 |
| ρ [g/cc] (caja de referencia) | 7.471 (caja promedio de 50 ps de NPT) | 7.554 | 7.554 |

- Isotropía excelente: (C11−C12)/2 = 19.39 vs C44 = 19.39 a 300 K y 23.02 vs 23.14 a 0 K; índice A_U ≈ 0.0006–0.0008; Voigt = Reuss = Hill a 3 cifras.
- Errores: a 300 K el error estándar por bloques de 5 ps da ±0.05–0.07 GPa por Cij (±0.02 GPa en K y G, ±0.1 K en ΘD), **pero subestima**: la asimetría |Cij − Cji| llega a 0.31 GPa en los términos cruzados (que deberían ser 0). Un error más realista para K y G es ≈ ±0.1–0.2 GPa, ≈ ±1 K en ΘD.
  A 0 K la diferencia entre δ = 1e-4 y 5e-4 es ≤ 0.2 GPa en cualquier Cij (≤ 0.4 % en E y G): resultado independiente de δ. Solo UNA muestra (un vidrio, una trayectoria): no hay barras entre muestras.
- Una de las 24 minimizaciones (δ = 1e-4, dir 2, signo −) paró por el tope de evaluaciones de fuerza (fmax 3.3e-4 eV/Å) y no por línea de búsqueda; sus presiones son antisimétricas con las del signo + a 0.1 bar, así que no sesga la columna. Las demás pararon por "linesearch alpha is zero".
- **Contraste con el piloto de tracción:** E = 54.8 GPa (Cij, 300 K) vs 54.6 GPa (ajuste de la curva a ε ≤ 0.01) y ν = 0.413 vs ≈ 0.41. La rigidez baja del piloto NO es un artefacto del ajuste ni de la tasa: es una propiedad de este vidrio con este potencial.
- Ablandamiento térmico (0 K → 300 K): E −15 %, G −16 %, K −6 %; el volumen de 300 K es 1.1 % mayor.
- **ΘD vs Méndez (312 K):** las constantes elásticas de este potencial dan ΘD = 209 K (300 K) y 228 K (0 K), mucho menos que 312 K. Para 312 K haría falta v_m ≈ 1.5 veces mayor (G del orden de 40 GPa). Con el mismo cálculo, módulos típicos de un vidrio Zr-Cu experimental (G ≈ 32, K ≈ 110 GPa, ρ ≈ 6.8 g/cc; valores de memoria, sin verificar) darían ≈ 270 K.
  Falta aclarar cómo se obtuvo el 312 K de Méndez (temperatura de Debye elástica, de calor específico o de ajuste a la densidad de estados de fonones): no son la misma magnitud. La ΘD de Anderson de acá es la elástica.

## Matriz de ensayos (2026-10-08)
Ver `matrix/README.md`: `in.tension` generalizado (eje x|y|z, rampa de calentamiento), 55 casos en `matrix/manifest.csv`, scripts de lanzamiento, estado y análisis en `matrix/scripts/`. Estado: pruebas cortas OK; ninguna etapa lanzada todavía.

## Cij por tasa de enfriamiento (Hill; 0 K con δ = 1e-4) — 2026-10-10
| muestra | K, G, E 300 K [GPa] | ν 300 K | ΘD 300 K [K] | E por tracción (media x,y,z, ε ≤ 0.01) | K, G, E 0 K [GPa] | ΘD 0 K [K] |
|---|---|---|---|---|---|---|
| 2e10 | 105.3, 19.39, 54.8 | 0.413 | 209.3 | 54.3 | 112.0, 23.09, 64.8 | 227.7 |
| 1e11 | 104.9, 18.35, 52.0 | 0.417 | 203.8 | 51.9 | 111.7, 22.47, 63.2 | 224.8 |
| 1e12 | 104.6, 17.49, 49.7 | 0.421 | 199.1 | 50.2 | 111.2, 21.77, 61.3 | 221.3 |
Isotropía (C11−C12)/2 ≈ C44 en todos; E de Cij y E de tracción coinciden a ≤ 1 % a 300 K. G baja ≈ 10 % de 2e10 a 1e12 (K casi no cambia). Ningún valor se acerca a los 312 K de Méndez.
A 0 K, δ = 5e-4 vs 1e-4: G cambia 0.4 % (2e10), 0.4 % (1e11) y 1.4 % (1e12).
