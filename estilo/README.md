# estilo/ — criterio único de figuras

- `bmg_estilo.py`: paleta, tamaños, rcParams, etiquetas, coma decimal, `guardar()` con revisión
  automática e índice de figuras.
- `muestra_estilo.py`: genera `/scratch/negishi/aaabt/bitacora/figs/estilo_{muestras,especies,cmapT}.png` para ver
  el estilo. Volver a correrlo después de tocar `bmg_estilo.py`.

Las reglas (y el porqué de cada una) están en la skill
`/scratch/negishi/aaabt/.claude/skills/figuras-bmg/SKILL.md` (fuera del repo). Resumen:

| Decisión | Valor | Por qué |
|---|---|---|
| Fuente | Times New Roman → Liberation Serif en Negishi; fórmulas STIX | es la fuente del informe de avance (usuario, 2026-11-13) |
| Decimales | coma, sin separador de miles | pedido del usuario; `locale` es_AR daría `1.000` |
| Paneles | uno por archivo | los paneles se arman en LaTeX (usuario, 2026-11-13) |
| Colores | Okabe-Ito, uno fijo por muestra | antes 1e12 era rojo en Tg y verde en mecánica |
| Referencias externas | negro/gris (`mendez`, `ryltsev`, `experimento`; conteos bibliográficos `biblio_mg` negro, `biblio_bmg` gris, línea sin marcador) | que no compitan con los colores de las muestras propias |
| Tamaño | `media` 7,2×5,4 cm (default) | incluir a tamaño natural en Beamer 16:9 |
| Salida | solo PNG 300 dpi en `../../bitacora/figs/` + `indice_figuras.csv` | sin PDF por figura (usuario); trazabilidad figura → script → datos |

Python: `/scratch/negishi/aaabt/venv_structural/bin/python` (matplotlib 3.11).
