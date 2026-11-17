# Fase 5: figuras, tablas y métodos de transporte

Todo se regenera desde los resultados ya calculados en `transporte/` (no corre simulaciones).

| archivo | qué es | se genera con |
|---|---|---|
| `fig_kappa_T.pdf/png` | (a) κ(T) para las 3 tasas, Tg de cada tasa, ajuste de recta quebrada y NEMD; (b) meseta, máximo de κ(t) y Méndez et al. | `figuras.py` |
| `fig_eta_D_T.pdf/png` | (a) η(T) con VFT y Arrhenius; (b) D por especie; (c) Stokes-Einstein | `figuras.py` |
| `tabla_kappa.tex` | κ(T) por tasa | `tablas.py` |
| `tabla_eta_D.tex` | η por tasa y D por especie en el líquido | `tablas.py` |
| `tabla_corridas.csv/.tex` | parámetros de cada corrida, leídos de los `run_*.sbatch` y salidas de cada `cal_*` (para apéndice) | `tablas.py` |
| `metodos_transporte.tex` | sección de métodos y limitaciones (fragmento para `\input`) | a mano |
| `vista_previa_metodos.tex/pdf` | documento mínimo que compila el fragmento con las tablas y figuras (8 páginas) | `pdflatex` (módulo `texlive/20220321`) |

Uso, desde `transporte/fase5/`:

```bash
/scratch/negishi/aaabt/venv_structural/bin/python figuras.py
/scratch/negishi/aaabt/venv_structural/bin/python tablas.py
```

Fuentes que leen los scripts: `../kappa_tasas.csv`, `../eta_T.csv`, `../eta_tasas.csv`, `../eta_largo.log`, `../eta_vft.txt` y `../../structural_analysis/Tg_resultados.csv`. Escritos a mano en `figuras.py`: κ de NEMD (`../nemd/README.md`) y la Tabla 3 de Méndez et al.

## Pendientes del texto

- Citas: las claves están listadas al principio de `metodos_transporte.tex`. Las marcadas VERIFICAR hay que confirmarlas antes de compilar con el `.bib`, en particular `sheng_eam`, `surblys`/`boone` (si aplican al EAM) y `zhang_maginn`.
- Completar cómo se generó la configuración fundida (`in.melt` de la referencia a 1e12) y la duración del equilibrado a 2000 K (`restart.equil2000K`).
- `\ref{sec:tg}` apunta a la sección de Tg de la tesis: ajustar la etiqueta.
- La semilla de cada réplica sale de la convención de los lanzadores; no queda en los logs.
