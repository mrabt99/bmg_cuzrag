# Etapa 4 — RDF y BADF (réplica de Figs 2, 3a y 4 de Ryltsev et al., JCP 145, 034506, 2016)

Cu-Zr-Ag N=65536, 3 tasas (2e10, 1e11, 1e12 K/s), 17 T por tasa. Ag es especie propia en todo
(en el RDF total y en los vecinos del BADF). No se tocó `boo/` ni `aggregate/`.

## Mapeo config → T
`NN = 01..17`, `T = 1900 − 100·(NN−1)` K, `NN=17` = vidrio 300 K (igual que Etapas 2/3).
En 1e12 no hay restart de sweep: `01` = seg1, que termina a 1900 K.

## Manifests corregidos
Los `<tasa>/manifest.txt` originales apuntan a `/scratch/negishi/aaabt/65536_*/`, que ya no existe
(reorg del 2026-09-18). `manifests/<tasa>.txt` son copias con las rutas nuevas
(`65536_1e10 → bmg_cuzrag/rate_2e10`, `65536_1e11 → rate_1e11`, `65536_1e12 → rate_1e12_reference`);
los 51 archivos existen. Los originales no se modificaron.

## Tarea 1 — RDF (`rdf/`)
- `run_rdf.lmp` + `batch_rdf.sbatch` (standby, 8 cores): `read_restart` desde el manifest corregido,
  `pair_style eam/alloy` con `rate_2e10/ZrCuAg.lammps.eam` (orden Zr Cu Ag), `run 0`.
- **r-grid (idéntico en las 3 tasas y las 17 T):** `compute rdf 500 cutoff 10.0` → 500 bins en 0–10 Å,
  Δr = 0.02 Å, r = centro de bin (0.01 … 9.99 Å). Total all-all. Requiere `comm_modify cutoff 12.5`
  (cutoff del rdf + skin 2.0).
- `Ncum` = coordinación acumulada N(<r) por átomo, tal cual la entrega LAMMPS.
- `rdf_all.csv`: rate, T, r, g_total, Ncum (25 500 filas). `rdf_qc.csv`: caja L, ρ_num y
  N(<9.99 Å) vs ρ·(4/3)π·r³: cociente 0.9967–1.0002 en las 51 configs. La caja es cúbica; ρ_num
  0.0520–0.0570 Å⁻³ (1900→300 K), consistente entre tasas.
- `build_rdf_csv_and_figs.py` genera los CSV, `fig_A_rdf_heatmap.png` y `fig_B_rdf300_overlay.png`.
  Fig A usa celdas discretas de ±50 K (17 T, sin interpolar). Fig B: eje derecho N(<r) con ticks 12 y 24.
- El eje r está en Å (Ryltsev usa unidades reducidas, pico ≈ 0.06 ↔ ≈ 2.8 Å aquí).

## Tarea 2 — BADF (`badf/`)
- Configs: `../<tasa>/snaps/17.dump` (vidrio 300 K) y `08.dump` (**líquido de referencia = 1200 K en
  las 3 tasas**, ≈ 500 K sobre Tg ≈ 660–690 K). Se verifica que la caja del dump coincida con `rdf_qc.csv`.
- **Cutoff = primer mínimo de g_total(r) de esa config** (de `rdf_all.csv`): Savitzky-Golay
  (9 bins = 0.18 Å, polinomio 3); pico = máximo del suavizado en 2.0–3.6 Å; mínimo = primer mínimo
  local a r > r_pico + 0.2 Å. Se usa el g(r) total y se aplica a todas las especies (también para
  centros Cu). Valores (Å): ver `badf_meta.csv`; vidrio 3.75 (1e12) / 3.79 (1e11) / 3.73 (2e10);
  líquido 3.93 / 3.83 / 3.89. El mínimo del vidrio es plano (±0.1 Å) y el cutoff cambia <CN> ≈ 0.2.
- Ángulos: `cKDTree` periódico, vectores en mínima imagen, todos los pares de vecinos de cada átomo
  central. `center=all` (65 536 centros) y `center=Cu` (type 2, 31 130 centros).
- **Normalización:** histograma de 1° (180 bins en 0–180°, `alpha` = centro de bin) de todos los
  ángulos acumulados; `P(α_k) = 180·cuentas_k / N_total`. Es una densidad respecto de x = α/π con
  ∫₀¹ P dx = 1 (media de P sobre 0–180° = 1; distribución uniforme → P = 1). Es la escala de la
  Fig 4 de Ryltsev. Por radián: `P/π`. Sin peso 1/sin α. Misma función para vidrio, líquido y refs.
- Referencias ideales (`badf_refs.csv`): icosaedro (12 vértices, permutaciones de (0,±1,±φ)) y hcp
  (c/a = √(8/3), 12 vecinos), mismo histograma. Al ser ideales son deltas en bins de 1° (picos
  63.4°/116.6°/180° para ico); **en `fig_C_badf.png` se muestran con un ensanchamiento gaussiano
  σ = 3° solo visual** (el CSV es el histograma exacto).
- Fig C: dos paneles (todo átomo central / Cu-centrado). Líquido 1200 K en gris (las 3 tasas casi
  superpuestas).

## Observaciones
- A 300 K las 3 tasas difieren poco: g_max/g_min = 8.06 / 8.47 / 8.67 (1e12 / 1e11 / 2e10),
  N(<r_min) = 12.8 / 13.0 / 13.2, tendencia monótona con la tasa. Nuestra tasa más lenta (2e10) es
  ~13× más rápida que la más lenta de Ryltsev (1.5e9), donde ellos ven la diferencia clara.

## Tarea 3 — RDF parcial por especie (`rdf_partial/`)
Pedido de Jaroszewicz (2026-09-23) tras ver la presentación: probar el RDF por especie, esperando
ver la separación del 2do pico que él y Mendez reportan en su Figura 5
(`../../referencias/Fugura5_Mendez.png`, heatmap g_total(r) vs T con el 2do pico partiéndose en dos
franjas cerca de Tg) — pero centrando en Cu y en Zr.

- Mismo esquema que Tarea 1 (`run_rdf_partial.lmp`, mismos 51 restarts vía `../manifests/`), pero
  `compute rdf 500 * * 1 * 2 * cutoff 10.0`: 3 pares tipo↔tipo, `*` = wildcard "todos los tipos"
  (no `0`; la sintaxis vieja con `0` no existe en esta versión). Da 3 g(r): total (all-all, igual a
  Tarea 1), centrado en Zr (type 1, vecino cualquiera), centrado en Cu (type 2, vecino cualquiera).
  Composición fija en las 51 configs: nZr=31129, nCu=31130, nAg=3277.
- `build_rdf_partial_csv_and_figs.py` → `rdf_partial_all.csv` (25 500 filas) + `fig_D_rdf_heatmap_Zr.png`
  + `fig_E_rdf_heatmap_Cu.png` + `fig_F_rdf300_species_zoom.png`. Sanity: g_total coincide con
  `rdf/rdf_all.csv` (max|diff| ~5e-7, redondeo de escritura del CSV).
- **Resultado: esta descomposición (centro=especie, vecino=cualquiera) NO muestra la separación del
  2do pico** — tanto centrado en Zr como en Cu, el 2do pico se ve como una sola joroba ancha en las
  17 T y las 3 tasas, igual que en el g_total de Tarea 1. Posibles causas a explorar: (a) puede que
  la Fig 5 de Mendez use la descomposición por *pares del mismo tipo* (g_Cu-Cu, g_Zr-Zr) en vez de
  especie-centro-vs-todos; (b) nuestro g(r) es de un único snapshot por T (`run 0`, sin promediar en
  el tiempo) — más ruido que podría estar enmascarando una separación sutil que ellos sí resuelven
  promediando varios frames.

## Tarea 4 — RDF de pares mismo tipo, Zr-Zr y Cu-Cu (`rdf_sametype/`)
Segundo intento tras el resultado negativo de la Tarea 3, probando la hipótesis (a): descomponer
por pares del mismo tipo en vez de especie-centro-vs-todos.

- Mismo esquema (`run_rdf_sametype.lmp`, mismos 51 restarts): `compute rdf 500 * * 1 1 2 2 cutoff 10.0`
  → total (sanity, coincide con Tarea 1), Zr-Zr, Cu-Cu.
- `build_rdf_sametype_csv_and_figs.py` → `rdf_sametype_all.csv` (25 500 filas), `fig_G_rdf_heatmap_ZrZr.png`,
  `fig_H_rdf_heatmap_CuCu.png`, `fig_I_rdf300_sametype_zoom.png`.
- **Resultado: POSITIVO para Cu-Cu.** A 300 K, g_Cu-Cu(r) tiene dos sub-picos bien definidos en
  r≈4.4-4.5 Å y r≈5.4 Å con un valle marcado en r≈5.0 Å (`fig_I_rdf300_sametype_zoom.png`,
  panel derecho). g_Zr-Zr(r) sigue siendo una sola joroba ancha (sin split), igual que el total y
  que Tarea 3.
- Cuantificado en `rdf_sametype/CuCu_split_summary.csv` (peak1/dip/peak2/profundidad, 17 T × 3 tasas):
  la profundidad del valle (promedio de los dos sub-picos menos el valle) es chica y casi plana a
  T≥900-1100 K (~0.08-0.12), y crece fuerte por debajo de esa T hasta ~0.29-0.37 a 300 K, en las 3
  tasas por igual. La posición del 2do sub-pico salta de r≈5.0-5.2 Å a r≈5.4 Å entre T=900 K y
  T=700 K — coincide con la zona de Tg≈660-690 K ya estimada (`compute_tg.py`).
  Es la misma "emergencia brusca cerca de Tg" que describe la Fig 5 de Mendez, pero aquí
  específicamente en la correlación Cu-Cu.
- `fig_J_CuCu_split_emergence.png` (tasa 2e10, curvas de g_Cu-Cu(r) coloreadas por T de 1900 K
  a 300 K): visualiza la transición de una joroba única (violeta, alta T) a dos picos separados
  (amarillo, 300 K).
- Conclusión provisoria: la separación del 2do pico SÍ aparece por especie, pero es específica del
  par Cu-Cu (no Zr-Zr, no el "centrado en Cu vs todos" de Tarea 3). Consistente con que Cu (radio
  menor, minoría) sea el que arma entornos icosaedrales/SRO más definidos en este BMG.

## Tarea 5 — RDF cruzado Zr-Cu (`rdf_cross/`)
Completa la descomposición binaria: falta el par cruzado además de Zr-Zr y Cu-Cu (Tarea 4).

- `compute rdf 500 * * 1 2 2 1 cutoff 10.0`: total (sanity) + Zr-centro/Cu-vecino + Cu-centro/Zr-vecino.
  Sanity extra: g_ZrCu(r) == g_CuZr(r) exactamente (max|diff|=0, como debe ser por definición de
  par cruzado) — confirma que la normalización del `compute rdf` de LAMMPS es simétrica.
- `build_rdf_cross_csv_and_figs.py` → `rdf_cross_all.csv`, `fig_K_rdf_heatmap_ZrCu.png`,
  `fig_L_rdf300_cross_zoom.png`, `fig_M_ZrCu_cross_emergence.png` (análoga a `fig_J` de Cu-Cu:
  g_Zr-Cu(r) coloreada por T, tasa 2e10).
- **Resultado: parcial/más débil que Cu-Cu, pero con la misma tendencia.** g_Zr-Cu(r) tiene un pico
  principal (r≈4.8-5.1 Å) y un hombro (r≈5.5-5.8 Å) separados por un valle mucho menos marcado que
  en Cu-Cu — a simple vista, con las ventanas de búsqueda ajustadas a esta forma (`ZrCu_split_summary.csv`,
  pico1 4.4-5.15 Å / valle 5.15-5.55 Å / pico2 5.5-6.15 Å; **ojo:** las ventanas de `CuCu_split_summary.csv`
  no sirven acá, la forma es distinta), la profundidad del valle crece de forma **gradual y monótona**
  con el enfriamiento: ~0.03-0.05 a 1900 K → ~0.09-0.11 a 900 K → ~0.28-0.29 a 300 K (las 3 tasas).
  A diferencia de Cu-Cu, que tiene un salto neto de posición del 2do sub-pico entre 900 K y 700 K
  (fingerprint más agudo cerca de Tg), acá el crecimiento es suave, sin discontinuidad clara en Tg.

## Conclusión de las Tareas 3-5 (descomposición por especie)
De las 4 descomposiciones probadas, **Cu-Cu** es la que muestra el fingerprint más nítido: split neto
del 2do pico con un salto brusco de posición justo cerca de Tg (900→700 K). **Zr-Cu cruzado** muestra
una versión más débil y gradual del mismo efecto (el valle se profundiza al enfriar pero sin el salto
discontinuo). **Zr-Zr** y el centrado-especie-vs-todos (Tarea 3) no muestran nada apreciable. Conclusión:
el efecto es más fuerte cuanto más involucrado está el Cu (Cu-Cu > Zr-Cu > Zr-Zr ≈ 0), consistente con
que Cu (radio menor, especie minoritaria) sea el que arma el orden icosaedral/SRO más definido en este
BMG.

## Tarea 6 — Figuras estilo Figura 4 de Mendez, por par (`fig4_style/`)
Pedido explícito: reproducir el estilo de `../../referencias/Figura4_Mendez.png` (panel izq.: g(r)+offset
apiladas por T, waterfall; panel der.: zoom del 2do pico sin offset) para Cu-Cu, Zr-Zr y Zr-Cu.

- `build_fig4_style.py` lee `rdf_sametype/rdf_sametype_all.csv` (Cu-Cu, Zr-Zr) y
  `rdf_cross/rdf_cross_all.csv` (Zr-Cu), tasa 2e10, las 17 T. Offset = 0.45/curva (panel izq.),
  colores con `nipy_spectral` (T baja = rojo, T alta = violeta, igual sentido que la referencia).
  Panel derecho usa colorbar en vez de etiquetas por curva (las 17 líneas convergen a la derecha y
  se superponían los textos).
- Salida: `fig_N_fig4style_CuCu.png`, `fig_O_fig4style_ZrZr.png`, `fig_P_fig4style_ZrCu.png`.
- Confirma visualmente lo ya cuantificado en Tareas 4-5: en Cu-Cu el panel derecho muestra el vidrio
  (rojo) con dos picos netamente separados (~4.4 y ~5.4 Å, valle en ~5.0) mientras el líquido a alta T
  (violeta) es una sola joroba; en Zr-Zr no hay separación en dos picos comparables (un pico principal
  con hombros menores, forma estable en T); en Zr-Cu se ve una versión intermedia -- un hombro que se
  define un poco más al enfriar pero de forma más gradual y menos marcada que en Cu-Cu.
