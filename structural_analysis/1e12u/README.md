# 1e12u: análisis estructural de la 1e12 UNIFICADA (`rate_1e12/`, 17/11/2026)

`1e12u` es la etiqueta de la 1e12 rehecha con el mismo protocolo que 2e10 y 1e11. Las carpetas `1e12/` y los CSV sin sufijo siguen siendo de la referencia vieja y **no se tocaron**.

| paso | qué se hizo | salida |
|---|---|---|
| thermo + Tg | `extract_thermo.py --style avetime`, `tg_1e12u.py` (compute_tg.py sin cambios) | `thermo_1e12u.csv`, `Tg_1e12u.csv` (Tg 698.0 / 698.0 K) |
| manifest | 17 restarts (sweep + seg1..15 + glass300K) | `manifest.txt` (copia en `etapa4/manifests/1e12u.txt`) |
| dumps | `../batch_restart_to_dump_pmi2.sbatch 1e12u` (job 46084063) | `snaps/NN.dump` |
| BOO | `../boo/batch_boo_pmi2.sbatch 1e12u <pot>` (job 46084064) | `../boo/1e12u/NN.boo.dump` |
| Voronoi | `voronoi_1e12u.py` (OVITO default, `module load anaconda/2024.02-py311`) | `abundance_voronoi_1e12u.csv`, `cn_distribution_voronoi_1e12u.csv` |
| tablas | `tabla_1e12u.py` (`--validar` reproduce las razones viejas 2e10/1e12 exactas) | `tabla_compare_rates_u.csv`, `tabla_regimen_pool_u.csv` |
| RDF | `etapa4/rdf*/batch_*.sbatch 1e12u` (jobs 46084065-68), `etapa4/build_1e12u_csv.py` | `etapa4/rdf*/..._1e12u.csv` |
| BADF | `etapa4/badf/badf_1e12u.py` | `badf_1e12u.csv`, `badf_meta_1e12u.csv` |

## Resultado: razones 2e10/1e12u (vítreo 500-300 K)

| | 1e12u | 1e12 vieja |
|---|---|---|
| Cu-ico BOO | 1.61 ± 0.06 | 1.74 ± 0.07 |
| Cu-ico Voronoi | 1.28 ± 0.02 | 1.29 ± 0.02 |
| Zr Z14 / Z15 / Z16 | 1.12 / 1.20 / 1.22 | 1.14 / 1.21 / 1.22 |
| Ag-ico Voronoi | 1.20 ± 0.10 | 1.06 ± 0.08 |

En el líquido todo da compatible con 1. Desdoblamiento Cu-Cu a 300 K: 0.263 (antes 0.252). g_max/g_min a 300 K: 8.06 (antes 7.86). Las conclusiones no cambian. La bitácora (secciones 3-6) ya usa la 1e12u.
