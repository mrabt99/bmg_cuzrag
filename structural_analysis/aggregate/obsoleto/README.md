# Figuras obsoletas (movidas el 2026-09-21)

No usar en presentaciones ni tesis. Ver `../auditoria_Z_CN_vs_indice/INFORME.md`.

| Archivo | Por qué es obsoleto | Reemplazo vigente |
|---|---|---|
| `compare_n_Z14.png`, `_Z15.png`, `_Z16.png` | `Coordination == 14/15/16` (no índice exacto) + `relative_face_threshold=0.01`; no es Frank-Kasper. Sus títulos aún dicen "Frank-Kasper ... (Voronoi CN==k)" porque son los PNG originales sin regenerar | `../figs/compare_rates_Zr_Z14/15/16.png` (índice exacto, OVITO default) |
| `compare_n_ico_voronoi.png` | ico ⟨0,0,12,0⟩ con umbral de cara (hasta ~27 % menos en Cu que con OVITO default) | `../figs/compare_rates_Cu_ico_Voronoi.png`, `..._Ag_ico_Voronoi.png` |
| `compare_ico_ratio_slow_fast.png` | la parte Voronoi usa el ico con umbral; la parte BOO no está afectada | sin reemplazo directo; razones vigentes en `../tabla_regimen_pool.csv` |

Aviso: `../step5_compare_rates.py` sigue escribiendo `compare_ico_ratio_slow_fast.png` en `aggregate/` si se reejecuta.
