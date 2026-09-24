# Auditoría de las curvas Z14/Z15/Z16 ("Voronoi CN==14/15/16")

Fecha: 2026-09-21. Alcance: todo lo que en `structural_analysis/` produce, transforma o grafica fracciones Z14/Z15/Z16 y icosaédricas por Voronoi. No se modificó ningún archivo existente; todo lo nuevo está en esta carpeta.

## 1. Respuesta corta

Las curvas etiquetadas "Frank-Kasper Z14/Z15/Z16 (Voronoi CN==14/15/16)" (`aggregate/compare_n_Z14.png`, `compare_n_Z15.png`, `compare_n_Z16.png`) usan la **definición 2, número de coordinación total**, y además una teselación con umbrales de cara. **No son comparables con la tabla citable de Etapa 3** (índice exacto ⟨0,0,12,k⟩) ni miden Frank-Kasper.

Hay dos pipelines Voronoi paralelos en el proyecto y solo uno es el que se cita:

| | Pipeline A (17-sep) | Pipeline B (18-sep) |
|---|---|---|
| Script | `aggregate/step2_voronoi.py` → `step3` → `step4` → `step5` | `voronoi/voronoi_51.py` → `step7/8/9/11/12` |
| Definición Z14/15/16 | `Coordination == 14/15/16` (l. 66-68) | índice exacto ⟨0,0,12,2/3/4⟩ (l. 34-35, 73-77) |
| Teselación | `edge_threshold=0.1`, `relative_face_threshold=0.01`; radios sin efecto (bug) | OVITO default |
| Salidas | `step2_voronoi.csv`, `tabla_clusters_51.csv`, `abundancia_*.csv`, `compare_n_*.png`, `compare_ico_ratio_slow_fast.png` | `voronoi/abundance_voronoi.csv`, `tabla_compare_rates_3tasas.csv`, `tabla_regimen_pool.csv`, `tabla_fk_zr.csv`, `resumen_etapa3_numeros.md`, `figs/compare_rates_*.png` |
| Estado | superado, pero sin marca de obsoleto | vigente |

## 2. Hallazgos

### H1. CN==k y el índice exacto miden poblaciones distintas (magnitud ~10-1000×)

Misma teselación (OVITO default), Zr, 2e10 K/s, 300 K:

| | índice exacto | CN == 12+k | átomos con CN==12+k que sí son ⟨0,0,12,k⟩ |
|---|---|---|---|
| Z14 | 1.63 % | 15.0 % | 10.8 % |
| Z15 | 2.80 % | 36.3 % | 7.7 % |
| Z16 | 1.45 % | 32.9 % | 4.4 % |

Sobre las 51 configuraciones (Zr), la mediana de la fracción de átomos con CN==14/15/16 que tienen el índice exacto es 3.1 % / 1.8 % / 1.1 % (máximos 10.8 / 7.8 / 4.5 %). El cociente CN/exacto en Zr tiene medianas 32 / 56 / 88. En Cu y Ag es de 2-3 órdenes de magnitud, porque el índice exacto casi no existe ahí.

La figura `fig_zr_definiciones_2e10.png` lo muestra: en Zr Z14, la curva exacta va de 0.13 a 1.6 % al enfriar; la curva de la figura actual va de 28 a 37 %.

Los valores exactos coinciden con los del piloto del proyecto (Zr 2e10 300 K: 1.63 / 2.80 / 1.45 %, chequeo de `step7`). Que ese piloto coincida con el paper no lo verifiqué; lo que sí está validado contra Mendez & Jaroszewicz (Fig. 6) es el CN medio por especie con la teselación default.

### H2. La razón 2e10/1e12 cambia de conclusión según la definición

Régimen vítreo (500/400/300 K, N≥20 por T), Zr:

| | índice exacto (tabla Etapa 3) | CN default | CN + umbral (figura actual) |
|---|---|---|---|
| Z14 | 1.141 ± 0.045 | 0.985 ± 0.012 | 0.991 ± 0.008 |
| Z15 | 1.211 ± 0.035 | 1.018 ± 0.008 | 1.046 ± 0.008 |
| Z16 | 1.223 ± 0.049 | 1.016 ± 0.008 | 1.066 ± 0.015 |

Con el índice exacto hay más orden a tasa lenta (razón >1 a 3-6σ). Con CN las curvas de las tres tasas se superponen (razón ≈ 1) y en la figura actual parecen mostrar "no hay dependencia con la tasa". Eso es una propiedad de la definición: CN==14 está dominado por poliedros genéricos que no distinguen el vidrio más o menos relajado. No hay que leerlo como resultado físico. Los errores de CN son minúsculos porque N es del orden de 10⁴, pero son errores Poisson (ver limitaciones).

### H3. Las curvas CN también dependen fuertemente de los parámetros de teselación

`relative_face_threshold=0.01` baja el CN medio de 13.99 a 12.89 y reasigna átomos entre clases de CN. Zr 2e10 300 K, N con CN==14 / 15 / 16:

- default: 4682 / 11305 / 10254
- con umbral: 11435 / 10606 / 3336

CN==16 pasa de 33 % a 11 % solo por el umbral. El índice exacto es mucho más estable (Zr Z14 506 → 423, Z16 451 → 254) pero tampoco es invariante. Esto significa que las curvas CN no solo miden otra cosa: su valor numérico depende de un parámetro de teselación.

### H4. Bug de radios en el pipeline A: confirmado empíricamente

`step2_voronoi.py` fija los radios de Slater sobre `pipeline.compute(0)` (una copia) y activa `use_radii=True`. Reproduje esto en 5 configuraciones (2 tasas, 3 temperaturas) y el resultado con "radios" es **idéntico bit a bit** al de solo umbrales, en todas las especies. Los radios nunca se aplicaron. La memoria del proyecto ya lo registraba; aquí queda comprobado con los dumps. `voronoi_51.py --radii` (que sí fija los radios en `pipeline.source`) es la versión correcta para esa variante.

### H5. Señales de que estas figuras siguen "vivas" sin advertencia

- `step5_compare_rates.py` l. 77-80: el comentario dice que las figuras `n_ico_voronoi`, Z14, Z15, Z16 "ya estan bien" y las excluye del loop. Reejecutar step5 no las regenera ni las corrige, y el comentario da a entender lo contrario de lo que se encontró.
- `step5_compare_rates.py` l. 30-32: los títulos dicen "Frank-Kasper Z14 (Voronoi CN==14)". El título es honesto respecto de CN, pero "Frank-Kasper" no lo es.
- `step3_merge_validate.py` valida "Zr domina Z15 y Z16" y "Ag intermedio en Z14". Con CN eso se cumple casi por construcción (Zr tiene el CN medio más alto), así que no valida que la definición sea Frank-Kasper.
- Ningún README ni nombre de archivo marca `compare_n_Z*` como obsoletos.

### H6. Figuras ico Voronoi del pipeline A también están desactualizadas

`compare_n_ico_voronoi.png` y la parte Voronoi de `compare_ico_ratio_slow_fast.png` usan `N_ico_Voronoi` con umbral. Frente al default (`abundance_voronoi.csv`), la fracción con umbral / default es, sobre las 51 configs: Cu 0.73-0.97 (mediana 0.88), Zr 0.74-1.38, Ag 0.55-1.67 (Ag tiene N pequeño). Es un efecto de hasta ~27 % en Cu. La versión vigente es `figs/compare_rates_Cu_ico_Voronoi.png` y `..._Ag_ico_Voronoi.png`. La parte BOO no depende de Voronoi y no se ve afectada.

## 3. Lo que sí está bien (verificado)

- **Cadena de tablas**: `step2_voronoi.csv` = `tabla_clusters_51.csv` (153/153 filas, 0 diferencias en N_Z14/15/16). `abundancia_*.csv` = N / N_total_especie (0 discrepancias). Es decir, la figura CN grafica fielmente lo que calculó step2; el problema es la definición, no un error de transcripción. Un punto concreto: el 36.7 % de Zr Z14 a 300 K de la figura equivale a 11435/31129 y lo reproduje con OVITO desde el dump.
- **Pipeline B reproducible desde cero**: recalculé las 51 configuraciones y `abundance_voronoi.csv` coincide con diferencia máxima 5×10⁻⁷ (redondeo a 6 decimales).
- **Criterio exacto**: `voronoi_51.py` compara solo las columnas 2..5 (n3..n6) y no exige n≥7 = 0. Comprobé que exigirlo no cambia ningún conteo en las 612 celdas. Los índices ⟨0,0,12,2/3/4⟩ son los Z14/Z15/Z16 estándar de Frank-Kasper.
- **Consumidores del pipeline B**: step7, step8, step9, step11 y `etapa4/ovito/ovito_snapshot.py` declaran no usar CN== y sus entradas son `abundance_voronoi.csv`. Leí los scripts y sus fuentes de datos; no reejecuté sus cálculos de razones salvo los de Zr Z14/15/16 (H2), que coinciden con `tabla_fk_zr.csv` (1.141 / 1.211 / 1.223).

## 4. Inventario: qué usar en la presentación

| Artefacto | Definición | ¿Usar? |
|---|---|---|
| `figs/compare_rates_Zr_Z14.png`, `_Z15.png`, `_Z16.png` | índice exacto, default | **Sí.** Z16 con la nota de degeneración (Lu et al. PRB 102, 184201) |
| `figs/compare_rates_Cu_ico_Voronoi.png`, `Ag_ico_Voronoi.png`, `Cu_ico_BOO.png` | exacto / BOO | Sí |
| `tabla_regimen_pool.csv`, `tabla_fk_zr.csv`, `resumen_etapa3_numeros.md` | exacto, default | Sí |
| `compare_n_ico_boo.png` | BOO | Sí (independiente de Voronoi) |
| `compare_n_Z14.png`, `_Z15.png`, `_Z16.png` | CN + umbral | **No** como Frank-Kasper |
| `compare_n_ico_voronoi.png`, `compare_ico_ratio_slow_fast.png` | ico con umbral | No; reemplazar por `figs/` |
| `step2_voronoi.csv`, `tabla_clusters_51.csv`, `abundancia_*.csv` | CN + umbral | No citar para Z ni para ico Voronoi |

Limitación del reemplazo: el índice exacto solo tiene estadística para **Zr** (Cu y Ag: N<20 por T casi en todo el rango; Ag pool Z14 41/29, Z15 16/12, Z16 1/1). No existe hoy un equivalente exacto de los paneles Cu/Ag de `compare_n_Z*`.

## 5. Recomendaciones

1. Sacar `compare_n_Z14/15/16.png` de la presentación y usar los tres `figs/compare_rates_Zr_Z*.png`.
2. Si se quiere mostrar la sensibilidad a la definición como resultado (es un mensaje válido), usar `fig_zr_definiciones_2e10.png` de esta carpeta.
3. Si aun así se quieren curvas por CN, rotularlas "CN = 14" (sin "Frank-Kasper"), con teselación default y explicando que CN no identifica el poliedro.
4. **Hecho (2026-09-21, con tu confirmación):** `compare_n_Z14/15/16.png`, `compare_n_ico_voronoi.png` y `compare_ico_ratio_slow_fast.png` movidos a `aggregate/obsoleto/`, con un README allí. Los PNG conservan sus títulos originales ("Frank-Kasper ... CN==k"); no se regeneraron.
5. **Hecho:** en `step5_compare_rates.py` se corrigió el comentario de l. 77-80 (ahora dice que son CN con umbral, obsoletas, no Frank-Kasper) y los títulos `CLUSTER_TITLE` de Z14/15/16 (ahora "Coordinacion CN==k (Voronoi, umbral de cara; no es Frank-Kasper)"). Solo cambian texto y comentarios; la lógica no se tocó. Pendiente conocido: si se reejecuta step5, vuelve a escribir `compare_ico_ratio_slow_fast.png` en `aggregate/` (se dejó un comentario de advertencia allí).
6. Opcional: regenerar los paneles Cu/Ag exactos con estadística agregada (pool de T) en vez de por T, si se quiere mostrar algo para esas especies, con las salvedades de N pequeño.

## 6. Limitaciones y lo que NO se verificó

- No verifiqué la cita de Lu et al. sobre la degeneración de ⟨0,0,12,4⟩; se toma de la documentación del proyecto. Z16 exacto es un límite superior.
- No verifiqué el mapeo tipo→especie (1=Zr, 2=Cu, 3=Ag) de forma independiente; es coherente con que el CN medio por especie (Zr 15.4, Cu 12.6, Ag 14.0) reproduce Mendez & Jaroszewicz, pero no lo comprobé contra el `pair_coeff`.
- No auditamos BOO, ventanas de Tg, razones de Cu/Ag ico ni los valores de Etapa 4.
- Es un solo snapshot por (tasa, T), sin promedio temporal. Los errores son Poisson (√N), que subestiman por correlación espacial, como ya indica `resumen_etapa3_numeros.md`.
- No reejecuté `step5_compare_rates.py`, para no tocar sus PNG; verifiqué sus valores contra los datos.
- La comparación de la sección H1/H2 usa una sola configuración por T y tasa, igual que la tabla original.
- No encontré ninguna presentación (.pptx/.tex) dentro del proyecto; no sé qué figuras ya están en ella.

## 7. Contenido de la carpeta

| Archivo | Qué es |
|---|---|
| `INFORME.md` | este informe |
| `fig_zr_definiciones_2e10.png` | Zr Z14/15/16 vs T con las 3 definiciones (2e10) |
| `a1_control_ovito_5configs.py` | control OVITO: default vs umbral vs "radios" (bug) en 5 configs |
| `a2_recalculo_51configs.py` / `.csv` | recálculo de las 51 configs con OVITO default (exacto y CN) |
| `a3_integridad_y_razones.py` | integridad step2→tabla→abundancia y razones vítreas por definición |
| `a4_cn_vs_exacto_51configs.py` | comparación CN / exacto sobre las 51 configs |
| `a5_figura_evidencia.py` | genera la figura |

Ejecución: con `/apps/spack/negishi/apps/anaconda/2024.02-py311-gcc-8.5.0-lr57z2f/bin/python` (tiene OVITO 3.15.5 y pandas), desde esta carpeta. `venv_structural` no tiene pandas ni OVITO. Nota: a1 no guardó su salida; se reprodujo en consola.
