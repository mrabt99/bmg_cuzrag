# Etapa 3 — razón 2e10/1e12 por régimen (solo números)

Fuente: `tabla_regimen_pool.csv` (ico) y `tabla_fk_zr.csv` (FK exacto Zr). Scripts: `step8_*.py`, `step9_*.py`.
Voronoi = OVITO **default** (sin threshold), en toda la tabla. Ventanas compartidas 2e10/1e12:
vítreo 500/400/300 K, líquido 1900–900 K. Filtro N≥20 **por T** en ambas tasas; error = razón·√(1/N₁+1/N₂)
sobre el pool (piso Poisson, subestima por correlación espacial). Razón > 1 = más orden a tasa lenta.

## Régimen vítreo (500–300 K)

| Especie | Cantidad | Definición | Razón 2e10/1e12 | Pool 2e10/1e12 |
|---|---|---|---|---|
| Cu | ico | BOO (q6≥0.6, Ŵ6≤−0.16) | 1.745 ± 0.068 | 1797 / 1030 |
| Cu | ico | Voronoi ⟨0,0,12,0⟩ | 1.295 ± 0.021 | 8545 / 6599 |
| Zr | ico | Voronoi ⟨0,0,12,0⟩ | **no calculable** (2e10: 19/18/16 por T, <20) | — |
| Zr | FK Z14 | Voronoi ⟨0,0,12,2⟩ | 1.141 ± 0.045 | 1405 / 1231 |
| Zr | FK Z15 | Voronoi ⟨0,0,12,3⟩ | 1.211 ± 0.035 | 2592 / 2141 |
| Zr | FK Z16 | Voronoi ⟨0,0,12,4⟩ (límite superior*) | 1.223 ± 0.050 | 1357 / 1110 |
| Ag | ico | Voronoi ⟨0,0,12,0⟩ | 1.061 ± 0.084 | 332 / 313 |
| Ag | FK Z14/Z15/Z16 | Voronoi exacto | **no medible** (N<20 por T en los tres; pool Z14 41/29, Z15 16/12, Z16 1/1) | — |

## Régimen líquido (1900–900 K), ico

| Especie | Definición | Razón 2e10/1e12 | Pool | T retenidos | (razón−1)/error |
|---|---|---|---|---|---|
| Cu | BOO | 1.000 ± 0.236 | 36 / 36 | 1/11 (solo 900 K; ver nota 2) | 0.0 |
| Cu | Voronoi | 1.030 ± 0.023 | 4033 / 3917 | 11/11 | 1.3 |
| Zr | Voronoi | 0.934 ± 0.092 | 198 / 212 | 8/11 (descartados 1800, 1700, 1000) | −0.7 |
| Ag | Voronoi | 1.128 ± 0.167 | 97 / 86 | 3/11 (1200, 1000, 900) | 0.8 |

## De qué depende cada número

- La razón vítrea > 1 aparece en Cu-ico (BOO y Voronoi) y en Zr-FK (Z14/Z15/Z16); en los tres tipos de Zr-FK y en Cu-ico apunta en el mismo sentido (más orden a tasa lenta).
- Cu-ico vítreo depende de la definición: BOO 1.745 vs Voronoi 1.295.
- Zr-ico vítreo: sin dato (estadística, no física medida). Zr-FK: N de miles por tipo.
- Ag-ico vítreo: compatible con 1 (0.7σ). Ag-FK: sin estadística.
- Líquido: las cuatro razones están a ≤1.3σ de 1 (Cu-BOO: 1 solo T).

## Notas

1. *Z16 exacto ⟨0,0,12,4⟩ es degenerado en sistemas Cu-Zr (Lu et al., Phys. Rev. B 102, 184201, 2020): límite superior; pendiente cruzar con Z16 de BOO.
2. Cu-BOO líquido: el filtro N≥20 por T deja 1 de 11 T (900 K, 36/36). El resto tiene <20 en alguna tasa (0–3 átomos a 1900–1200 K, 2–15 a 1100–1000 K). La razón es de un solo T, en el borde frío de la ventana; valor previo sin filtro 1.218 ± 0.222 (67/55, ~5 por T). Lectura física: el orden icosaédrico se construye en el líquido sobreenfriado al acercarse a Tg; lejos de Tg se espera poca población (BOO: 0 a 1900–1700 K → 36 a 900 K), no es solo un problema de conteo. El criterio BOO es además más estricto que Voronoi ⟨0,0,12,0⟩ (107 átomos a 1900 K).
3. Ag-líquido: el filtro por T deja 3 de 11 T (todos ≤1200 K), no cubre la misma ventana que Cu.
4. Zr-líquido: alrededor del umbral (17–29 por T); 3 T descartados por 17, 19 y 18.
5. Filas Voronoi antes de esta corrección usaban `relative_face_threshold=0.01` + `edge_threshold=0.1` (respaldo: `tabla_regimen_pool_threshold001.csv`).
