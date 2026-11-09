# qgrad — pipeline de structural_analysis aplicado a la corrida graduada (Etapas 0–4 + búsqueda de Laves)

Corrida: `rate_gradiente_gammamin/` (12 segmentos). Perfil: 2000→1500 K = cadena 2e10 (compartida), 1500→700 K tasa
graduada 2e10→1.5e9 K/s (γmin de Ryltsev), 700→300 K a 1e11 K/s. **Es una prueba exploratoria.**
Todo lo nuevo vive acá; nada de `boo/ aggregate/ etapa4/ classify/ voronoi/ <tasa>/` se modificó (ver checksums, abajo).
Comparaciones = archivos NUEVOS con sufijo `_4curvas` (2e10, 1e11, 1e12, qgrad).

## Auditoría (Paso 0)
- 12 segmentos terminados limpios (sin fix halt, sin lost atoms/NaN; stderr = ruido PMIx/munge; único fallo real: intento pmix de
  seg02, ya reemplazado). Continuidad en bordes: ΔLx ≤ 0.013 Å, Δρ ≤ 0.002 g/cc, ΔEp ≤ 0.3 meV/at.
- `check_crystallization.py` por segmento (1..12) y global: sin señal conjunta caída-de-Ep + salto-de-ρ a 8σ.
- Mapeo T (NN = 1900−100(NN−1)): **01–05 (1900–1500 K) compartidos con 2e10** (mismo restart, no se recalculan; el 05 es
  sha256-idéntico a 2e10 seg4). 06–16 = seg01…seg11 (1400…400 K). 300 K: **17r** = vidrio relajado 10 ns NPT (`relax300K/`, job 43581066,
  8 h 47 min, limpio; equivalente del `glass300K` de las otras tasas) — **es el que usan todas las tablas y figuras**; 17u = seg12 sin
  relajar, se conserva como variante. T instantáneo al final de cada segmento: hasta 6 K del nominal (fluctuación de N=65536).
- Cada segmento guardó además 11 frames de dump (~cada 10 K, con velocidades) en `rate_gradiente_gammamin/*.dump.segNN`;
  NO se usaron (los restarts alcanzan; Laves da 0 a 300 K y un grano no se disuelve al enfriar).

## Qué se corrió (mismos criterios que las 3 tasas; scripts reusados o copiados con sufijo)
| Etapa | Qué | Dónde |
|---|---|---|
| 0 | restart→dump (`restart_to_dump.lmp` sin editar, wrapper `batch_qgrad.sbatch`, pmi2, standby) | `snaps/`, `manifest_*.txt` |
| 1 | `extract_thermo.py` + Tg de dos rectas (lógica y ventanas de `compute_tg.py`: liq 900–1200, vit 300–600, bines 10 K; reproduce exacto los Tg de las 3 tasas) | `thermo/` |
| 2 | BOO: `validate_boo.lmp` sin editar (EAM real, orientorder nnn 12/14/15/16, wl/hat). Ico BOO q6≥0.6 & Ŵ6≤−0.16. Voronoi OVITO **default**. FK exactos ⟨0,0,12,k⟩. Especies separadas; BOO y Voronoi separados | `boo/`, `aggregate/voronoi_*` |
| 3 | `tabla_clusters_qgrad.csv`, `tabla_compare_rates_4curvas.csv`, `tabla_regimen_qgrad.csv` (razón qgrad/ref), figuras `_4curvas` | `aggregate/` |
| 4 | RDF total/parcial/mismo-tipo/cruzado, heatmaps, overlay 300 K, BADF, Cu-Cu, estilo Fig.4 (copias parcheadas `_4curvas` de los builders; raw de las 3 tasas por symlink de solo lectura) | `etapa4/` |
| 5 | Laves | `laves/` |

**Criterios Etapa 3:** N≥20 por T en ambas curvas; error √(1/N₁+1/N₂) (piso Poisson); ventana por intersección exacta;
líquido T>Tg+100, vítreo T<Tg−100. **Tg del perfil qgrad = 631.9 K** (V 633.5, Ep 630.2; sensibilidad a ventanas ±20–40 K, en
`thermo/Tg_qgrad_y_sensibilidad.csv`). Adición propia: los T≥1500 K de qgrad (= 2e10) se excluyen de todos los pools (no independientes).

## Resultados (300 K = 17r relajado; entre corchetes lo que daba con 17u sin relajar)
- **Ep(300 K), vidrios relajados 10 ns, `run 0` sobre `restart.glass300K_*` (mismo potencial; `thermo/ep_run0_resultados.txt`)**:

  | | 2e10 | 1e11 | 1e12 | qgrad |
  |---|---|---|---|---|
  | Ep (eV/at) | −4.837286 | −4.834328 | −4.829831 | −4.837928 |
  | Ep − Ep(2e10) (meV/at) | 0 | +2.96 | +7.46 | **−0.64** |
  | [promedio últimos 5 ns, `log.dat`] | −4.83718 | −4.83415 | (sin log.dat) | −4.83783 |

  El `run 0` es un snapshot (coincide con el promedio de 5 ns a ≤ 0.2 meV/at). A 700 K qgrad − 2e10 = −3.3 meV: el tramo rápido 700→300 K borra casi toda la
  ventaja lenta. Relajar 10 ns casi no cambia el vidrio de qgrad (Ep −4.83781 → −4.83793).
- Líquido: qgrad = 2e10 (Cu-ico Voronoi 0.999 ± 0.020; Zr-FK ≤ 1.4σ de 1) — control correcto.
- Vidrio (500/400/300 K), razón qgrad/2e10: Cu-ico BOO **1.193 ± 0.038** (5.1σ) [1.184], Cu-ico Voronoi **1.105 ± 0.017** (6.4σ) [1.109];
  Zr Z14/Z15/Z16 = 1.034/1.018/1.016 (≤ 0.9σ); Ag-ico 1.045 ± 0.080 (0.6σ). Contra 1e11: Cu-ico BOO 1.545, Voronoi 1.241; contra 1e12: 2.082 / 1.431.
  Más icosaedros de Cu que 2e10, sin más FK de Zr.
- Cu-Cu (desdoblamiento del 2º pico): a 400–600 K qgrad sigue a 2e10 (profundidad del valle 0.339/0.284/0.272 vs 0.343/0.290/0.268). A 300 K, con **un snapshot**
  (run 0 del restart): qgrad 0.400, 2e10 0.353, 1e11 0.369, 1e12 0.322. **Con 10 frames** de cada relajación (`etapa4/rdf_sametype/cucu_multiframe.py`, g(r) por conteo
  de pares periódico, validado contra LAMMPS en 1e12 — el único cuyo último frame coincide con el restart: Δg ≤ 0.011, depth 0.3199 vs 0.3216):

  | | profundidad media de 10 frames (± std entre frames) | depth del g promediado (± bootstrap) |
  |---|---|---|
  | 1e12 | 0.321 ± 0.011 | 0.285 ± 0.007 |
  | 1e11 | 0.353 ± 0.009 | 0.318 ± 0.005 |
  | 2e10 | 0.369 ± 0.017 | 0.329 ± 0.006 |
  | qgrad | 0.376 ± 0.014 | 0.344 ± 0.007 |

  qgrad − 2e10 = +0.007 ± 0.007 (media de frames, 1.0σ) / +0.015 ± 0.009 (g promediado, 1.7σ): **el 0.047 del snapshot único era ruido térmico**
  (std por frame 0.014–0.017; la diferencia de dos frames sueltos tiene σ ≈ 0.02). Con multi-frame el orden es monótono en la lentitud
  (1e12 < 1e11 < 2e10 ≲ qgrad), qgrad compatible con 2e10 (+1–1.7σ, no concluyente). Los errores son SOLO ruido térmico entre frames; no incluyen la variación
  entre vidrios independientes (una sola muestra por tasa), así que subestiman. Posición del 2º sub-pico y salto entre 900 y 700 K: sin cambio. Ventanas de
  `CuCu_split_summary.csv` recuperadas por búsqueda exhaustiva (el script original no está en el árbol; reproducen las 51 filas a 5e-5).
- BADF: el mínimo del vidrio es plano y el cutoff propio de qgrad relajado (3.81 Å) es 0.08–0.10 Å mayor que el de 2e10 (3.73) y que el de 17u (3.71);
  con cutoff propio ⟨CN⟩ = 13.09 (qgrad) vs 12.77 (2e10), diferencia que es pura consecuencia del cutoff. Con cutoff **único** (3.770 Å, media de los 4):
  ⟨CN⟩ = 12.904 (1e12), 12.922 (1e11), 12.927 (2e10), 12.934 (qgrad): dispersión 0.03 (`etapa4/badf/badf_sharedcutoff.csv`).
  P(60°) Cu-centrado: 3.707/3.789/3.811/3.883 (1e12/1e11/2e10/qgrad).

## Hipótesis (NO es un resultado): "tasa efectiva" de qgrad según Ep
Cálculo en `thermo/tasa_efectiva.py` (→ `tasa_efectiva.csv`). **Supuestos: Ep(300 K, vidrio relajado) lineal en log10(tasa) y ln(N de cada cluster) lineal en log10(tasa).**
- Pendiente 2e10→1e11: 2.958 meV / log10(5) = **4.23 meV/at por década** (con Ep de promedios 5 ns: 3.03 meV → 4.34; la cuenta original). Linealidad: 1e11→1e12 = 4.50 y 2e10→1e12 = 4.39
  meV/década (dentro de ~6 %).
- Ep(qgrad) − Ep(2e10) = −0.642 meV/at → qgrad ≈ **0.15 décadas más lento que 2e10** (0.152 con la pendiente 2e10→1e11; 0.143 con 1e11→1e12; 0.146 con 2e10→1e12),
  tasa efectiva ≈ 1.4e10 K/s.
- Predicción de la razón qgrad/2e10 (vidrio 500/400/300 K) con esa tasa efectiva (interpolación log-lineal entre 2e10 y 1e11; ajuste a las 3 tasas entre paréntesis) vs medido:

  | cluster | predicha | medida | medida − predicha |
  |---|---|---|---|
  | Cu-ico BOO | 1.058 ± 0.008 (1.051) | 1.193 ± 0.038 | +0.135 (+3.5σ) |
  | Cu-ico Voronoi | 1.025 ± 0.004 (1.023) | 1.105 ± 0.017 | +0.080 (+4.7σ) |
  | Zr Z14 | 1.010 ± 0.008 | 1.033 ± 0.039 | +0.024 (+0.6σ) |
  | Zr Z15 | 1.018 ± 0.006 | 1.018 ± 0.028 | 0.000 (0.0σ) |
  | Zr Z16 | 1.024 ± 0.009 | 1.015 ± 0.039 | −0.009 (−0.2σ) |

  (Errores: Poisson de los pools, piso. La predicción es poco sensible a Δlog: 0.10/0.15/0.20 décadas → 1.038/1.057/1.077 para Cu-ico BOO. Los 5σ de "1.19" son
  respecto de 1; respecto de la razón predicha el exceso es 3.5σ (BOO) y 4.7σ (Voronoi).)
- **Lectura posible (hipótesis, no probada):** la energía del vidrio de qgrad equivale a una tasa efectiva ~0.15 décadas más lenta que 2e10, y los Zr-FK (Z14/Z15/Z16) son
  consistentes con esa tasa efectiva; pero el orden icosaédrico de Cu es mayor que el que predice esa tasa. Una interpretación compatible: los icosaedros de Cu formados en la
  ventana lenta sobreviven al tramo rápido 700→300 K mientras la energía "rejuvenece". Dato compatible (no prueba): razón qgrad/2e10 de Cu-ico a 700 K, justo antes del tramo
  rápido, = 1.32 (BOO) / 1.18 (Voronoi); a 300 K = 1.18 / 1.10, o sea el exceso se conserva parcialmente (Zr-Z16: 1.13 a 700 K → 1.07 a 300 K). Alternativas no descartadas: Ep no
  es lineal en log γ para qgrad (historia no monótona), o Ep y orden icosaédrico responden distinto a la historia térmica. Una sola muestra por tasa.

## Laves (Cu2Zr)
**Paper (Ryltsev et al., JCP 145, 034506; arXiv 1604.07024):** aleación **binaria Cu64.5Zr35.5**, N=5000, EAM de Mendelev;
"nano-sized crystalline grain of Cu2Zr … structure of Cu2Mg Laves phase" (= **C15**, MgCu2) a 300 K y γmin=1.5e9 K/s, identificado
**visualmente** (distribución espacial de centros ico/Z16 + inspección de snapshots, Fig. 11), sin análisis cristalográfico cuantitativo.
Nuestro sistema es Zr47.5Cu47.5Ag5 (≈1:1, Cu2Zr pide 2:1) y otro potencial: no es una réplica exacta.

**Método (umbrales fijados con las referencias, NO con qgrad; controles 2e10/1e11/1e12 en las 17 T):**
referencias C15 (24 at/celda) y C14 relajadas con el EAM real (0 K y 300 K NPT): E = −4.6563 / −4.6547 eV/at, a(C15) = 7.098 Å.
Detectores: **D** = subred de Zr (Cu/Ag borrados) es diamante cúbico/hexagonal (OVITO IdentifyDiamond; recall 100 % en C15/C14);
**T** = Zr con Voronoi ⟨0,0,12,4⟩, exactamente 4 vecinos Zr en tetraedro (≤10°, p99 de las referencias; recall 99.5 %);
**B** = Mahalanobis (q4,q6,Ŵ6) a la nube de referencia (χ² 99 %); clusters = componentes conexas del grafo de vecinos Voronoi.

**Veredicto: NO hay evidencia de fase de Laves en qgrad.**
- (a) Termodinámica: sin escalón exotérmico ni caída de V(T) en 1500–700 K (qgrad ≡ 2e10 hasta ~850 K; después la diferencia
  es el efecto de tasa, continuo). ✗ no coincide.
- (b) Zr-Z16 ≈ 2e10 (razón vítrea 1.016 ± 0.039; Fig. `aggregate/figs/compare_rates_Zr_Z16_4curvas.png`); Cu-ico sube (+11 %/+19 %) pero eso es orden icosaédrico, no Laves. ✗ (Z16 no coincide)
- (c) Diamante en subred Zr: **0 átomos en las 64 configuraciones** (2e10, 1e11, 1e12 todas las T = 51; qgrad 06–16 + 17u + 17r = 13). Zr Laves-like (T&B):
  2e10 = 0 en las 17 T; qgrad = 0 con 17r (con 17u había 1 átomo aislado a 500 y 300 K: desapareció al relajar, era una fluctuación). Zr-Z16 con exactamente 4 Zr vecinos: 9 de 484 (17r), 14 de 451 (2e10). ✗
- (d) Clusters Laves-like: ninguno (0 átomos Laves-like en qgrad 17r). Zr-Z16 + Cu-ico (laxo): cluster mayor = **51 at** (qgrad 17r, 300 K) [47 con 17u] vs 30 (2e10), 29 (1e11), 16 (1e12),
  1.4 % de esos poliedros, composición 10 % Zr / 90 % Cu / 0 % Ag (Laves pediría 33 % Zr): agregado de icosaedros de Cu, no un grano. ✗
- (e) g(r) parciales: diferencia qgrad−2e10 (rms 0.025–0.039) del mismo tamaño que 1e11−2e10 y 1e12−2e10; sin picos nuevos. S(k): sin
  picos de Bragg nuevos en (111)/(220)/(311)/(222) (poco sensible, ver límite); el máximo global de qgrad 17r (40.75) está en el primer pico del vidrio
  (|k| = 2.71 Å⁻¹, no en un Bragg) con máx/media = 13.6 vs 12.1–12.7 en los controles: fluctuación. ✗
- (f) Composición local: Ag excluido de los clusters (0 % dentro, ~2 % en el entorno vs 5 % nominal); sin segregación 2:1.
- **Límite de detección** (esferas de C15/C14 termalizado implantadas en el vidrio 2e10 a 300 K, sin relajar): 143 átomos (R=8 Å) → 50 Zr
  detectados por D (C15; 38 en C14); 279 átomos → 82. S(k) solo ve grano ≥ ~2000 átomos. O sea: **cualquier grano ≳ 1.5 nm sería visto;
  no hay ninguno**. Un grano de pocas decenas de átomos o muy defectuoso no está excluido.
- Dump para OVITO (con flags D, T, Z16, CuIco, BigCl): `laves/candidatos/qgrad_17r_flags.dump` (también `qgrad_17u_…` y `2e10_17_…`) (no hay candidato a grano; sirven para ver los agregados de icosaedros).

## Limitaciones
1. **qgrad no tiene tasa única**: su Tg (631.9 K) es la del perfil; el codo cae cerca del salto de tasa (700 K), así que mezcla historia lenta y rápida.
   Además, comparar "vidrio qgrad" vs "vidrio a tasa X" no es comparar tasas: la estructura final es el vidrio de ~Tf(1.5e9 K/s) + 400 K de enfriamiento rápido.
2. Puntos de T: los 12 coinciden nominalmente con la grilla; T instantáneo hasta 6 K de diferencia.
3. Los 5 puntos ≥1500 K son 2e10 (no independientes); el error Poisson es piso y subestima por correlación espacial (además comparten origen en 1500 K).
4. 300 K = 17r (relajado 10 ns), equivalente al `glass300K` de las 3 tasas originales. Los resultados con 17u se conservan en las tablas
   (config 17u) y en `features/qgrad/17u.npz`; las diferencias 17u→17r son menores (arriba).
5. Un solo snapshot por T (`run 0`) en RDF/BADF: más ruido que promediando frames (afecta sobre todo a la profundidad del valle Cu-Cu a 300 K, que ya se rehízo con 10 frames; ver arriba).
6. Paper binario Cu64.5Zr35.5 vs nuestro ternario ≈1:1; potencial distinto: ausencia de Laves acá no contradice al paper.
7. Voronoi de `tabla_clusters_qgrad.csv` es DEFAULT; `aggregate/tabla_clusters_51.csv` (original) usa radios+umbral: no comparar entre sí.
8. En `laves_analysis.py` no se hizo ajuste alguno con qgrad; T&B recupera ~97 % de átomos de referencia, no 100 %.

## Estado
Rehecho con 300 K relajado (17r) el 2026-09-30 con `run_17r.sh` (jobs 43610106–43610111). Único paso que falló en ese post-proceso: la exportación del dump de OVITO
(variable `$SA` sin definir en `post_17r.sbatch`); se hizo a mano y **`run_17r.sh` y `post_17r.sbatch` ya definen `SA`** (y se quitó el `|| true` que enmascaraba errores).
`laves_sq.py` se parchó para usar 17r. Añadidos después: Ep run 0 de las 4 curvas (job 43610965), Cu-Cu multi-frame, tasa efectiva.

## Aislamiento — checksums
`checksums.sbatch` (md5 de los 875 archivos de `structural_analysis/` excluyendo `qgrad/`): `checksums_antes.txt` (antes de empezar) vs
`checksums_despues.txt` (tomado al final de todo el trabajo, incluida la corrida con 17r, el Ep run 0 y el multi-frame; job 43611187):
```
$ diff checksums_antes.txt checksums_despues.txt        # (vacío: 875/875 idénticos; no apareció ningún archivo nuevo fuera de qgrad/)
```
Todo Python corrió con `PYTHONDONTWRITEBYTECODE=1`. Los scripts `run_tg_1e11.py` y los `step*` originales NO se ejecutaron (sobreescriben salidas);
se usaron copias `_qgrad`/`_4curvas`.
