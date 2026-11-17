# estructura_deformacion: relación entre la estructura inicial y los mecanismos de deformación

Carpeta creada el 17/11/2026 a pedido del usuario. Junta todo lo que asocia la estructura (Voronoi, BOO, icosaedros, Frank-Kasper) con la deformación (η de von Mises, D²min, STZ, nucleación). **Es un análisis en curso: no pasa a la bitácora hasta que esté maduro.**

## Qué hay y de dónde viene

- **Datos por caso:** se quedan en `mechanical/matrix/<caso>/analysis/` (`structure/`, `stz/`, `eta_loncha.npz`), al lado de los dumps. Los scripts los leen de ahí con la ruta absoluta `M`.
- `scripts/`: se mudaron con `git mv` desde `mechanical/matrix/scripts/`. Solo se cambiaron las rutas, no la lógica:
  - `M` = `mechanical/matrix` (casos, `manifest.csv`, `resumen.csv`);
  - `ED` = esta carpeta (`datos/`, `figs/`).

  Verificado: `nucleacion_figuras.py` reproduce idéntico `nucleacion_resumen.txt`.

| script | qué hace | salidas (`datos/`, `figs/`) |
|---|---|---|
| `structure_evolution.py` | Voronoi en los 41 cuadros de cada tracción (OVITO) | `<caso>/analysis/structure/voro_frames.npz` |
| `structure_vs_mechanics.py` | CN e icosaedros contra σ–ε | `estructura_vs_mecanica.csv`, `figs/estructura_*.png` |
| `stz_analysis.py` | D²min por ventana de Δε = 0.01, umbral p99.9 de la 1.ª ventana, clusters, control al azar | `<caso>/analysis/stz/stz_d2min.npz`, `stz_windows.csv` |
| `stz_figures.py`, `stz_largest_cluster_shape.py` | figuras de STZ y forma del cluster mayor | `stz_resumen.txt`, `stz_mayor_cluster.csv`, `figs/stz_*.png` |
| `nucleacion.py`, `nucleacion_figuras.py` | dónde nace la actividad: logística, riesgo relativo, celdas | `nucleacion_*.csv/txt`, `figs/nucleacion_resumen.png` |
| `evidencia_ico.py`, `ico_survival_eta.py` | pruebas A/B/C y supervivencia tipo Méndez | `evidencia_ico*.csv/txt`, `ico_survival_eta_1e9.csv` |
| `eta_loncha.py` (nuevo, OVITO) | η acumulada desde ε = 0 hasta 0.02 … 0.20, posiciones de ε = 0 | `<caso>/analysis/eta_loncha.npz` |
| `fig_mec_eta_loncha.py` (nuevo) | mapas de η a 0.06 y 0.10 y η por celdas contra fracción ico | `figs/mec_eta_*.png` |
| `evolucion_deformacion.py` (nuevo) | evolución con ε: mapas de STZ y de η sobre la loncha z < 12 Å (la de los mapas de q6), curvas contra ε | `datos/evolucion_deformacion.csv`, `figs/stz_loncha_*`, `figs/eta_loncha_*`, `figs/evol_*`, `figs/hoja_*` |

Uso: la mayoría corre con `/scratch/negishi/aaabt/venv_structural/bin/python`. Los que dicen OVITO necesitan `module load anaconda/2024.02-py311` y `python`.

## Primeros resultados (17/11, tracción en x, 300 K, 1e9 s⁻¹; 1e11, 2e10 y qgrad_ext; la 1e12u se sumó después, ver abajo)

- **La protección de los icosaedros depende de ε.** P(activo | Cu-ico) / P(activo | otro Cu), por ventana:
  - ≈ 0–0.4 antes del pico (ε ≤ 0.06);
  - ≈ 0.5–0.7 en el pico;
  - ≈ 0.8–1.1 en el flujo (ε ≥ 0.12).

  Los icosaedros frenan el *inicio* de la plasticidad. Una vez que el vidrio fluye, esa memoria se borra.
- η de los Cu-ico en clusters de 5 o más, frente a los otros Cu: mínimo ≈ 0.77–0.80 en ε ≈ 0.06–0.10; sube a 0.87–0.90 en 0.20.
- **Celdas de ~10.5 Å**, Spearman entre η y fracción ico: de −0.07 a −0.22. Lo más fuerte es qgrad_ext antes del pico.
- **STZ:**
  - activos: 0.2 % en ε = 0.03 → 2–3 % en 0.07 → 4–9 % en 0.10–0.20;
  - tamaño medio de cluster: 3–4 átomos al principio → 9–10 en el pico → 11–55 después;
  - qgrad_ext en ε = 0.10 tiene un cluster de 4759 átomos (avalancha que percola justo después del pico).

## Pendientes

- Plan completo (estadística, nucleación y coalescencia, fractura, 3D): `PLAN.md` (17/11/2026).

## 1e12u sumada (17/11/2026)
- Por caso (OVITO, log `datos/ovito_1e12u.log`): `structure_evolution.py --sample 1e12u` (15 casos: 300 K × 3 tasas, 450 y 540 K a 1e9), `stz_analysis.py` (mismos 15; filas agregadas a `datos/stz_windows.csv`, respaldo `stz_windows.csv.bak_pre_1e12u`), `eta_loncha.py 1e12u_300K_1e9_x`.
- Copias nuevas (originales intactos): `stz_figures_1e12u.py` → `datos/stz_resumen_1e12u.txt`, `figs/stz_*_1e12u.png` (reproduce idénticas las otras muestras); `stz_largest_cluster_shape_1e12u.py` (sin correr). `nucleacion.py --sample 1e12u --tag _1e12u` → `datos/nucleacion_casos_1e12u.csv`, `nucleacion_indices_1e12u.csv`.
- 1e12u (300 K, 1e9): activos estacionarios 3.49 % (vieja 3.62), tamaño medio 10.4 (10.6), mayor 149 (143); actividad relativa Cu-ico 0.58 (0.64), Cu no-ico 1.18, Zr CN15 0.82, Ag 1.12. AUC M4 0.593 (0.577); coef ico −0.13, entorno ico −0.16; riesgo relativo Cu ⟨0,0,12,0⟩ 0.42–0.56 según la tasa.
- **Los tres que faltaban, rehechos con 1e12u (17/11/2026, a pedido):**
  - `evidencia_ico_1e12u.py` (copia; solo procesa las muestras de SAMPLES) → `datos/evidencia_ico_1e12u.csv`, `evidencia_ico_resumen_1e12u.txt`, `figs/evidencia_ico_1e12u.png`. (B) pares emparejados 1e12u: razón 0.60 / 0.56 / 0.51 (1e8/1e9/1e10; vieja 0.51/0.47/0.46); las demás muestras 0.30–0.44. La protección local sigue en todos los vidrios y es la más débil en el rápido. (C) inicio 0.5 %: 1e12u 0.027 / 0.039 / 0.059 (vieja 0.026 / 0.034 / 0.061): sigue siendo el que arranca antes a 1e8 y 1e9, pero a 1e9 por menos (0.003 contra 1e11). (A, control qgrad exploratorio) igual que antes: Ep predice 2.320 (obs. 2.263), f_ico 2.439.
  - `structure_vs_mechanics_1e12u.py` (copia) → `datos/estructura_vs_mecanica_1e12u.csv`, `datos/estructura_resumen_1e12u.txt`, `figs/estructura_*_1e12u.png`. 1e12u: f0 = 0.0763 (vieja 0.0754), pérdida al pico 9.8 / 12.3 / 17.5 %, supervivencia a 0.20 ≈ 0.11. Correlaciones con 4 vidrios: σmax vs Ep r = −1.00, vs f_ico +0.99 (igual que antes).
  - `evolucion_deformacion.py` vuelto a correr tal cual (ya leía 1e12u); salidas anteriores en `respaldo_pre_1e12u_evolucion/`. 1e12u: P(activo|Cu-ico)/P(activo|otro Cu) 0.21–0.25 en ε 0.02–0.04, 0.48–0.49 en 0.05–0.06, 0.64–0.69 en 0.07–0.09 → ≈1 en el flujo (mismo patrón que 1e11); η de Cu-ico en clusters / Cu: mínimo 0.76 en ε = 0.06; Spearman de celdas −0.08 a −0.13. Único aviso de estilo: `evol_eta_celdas` (Spearman adimensional, benigno).
