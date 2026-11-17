# Plan: estructura inicial ↔ mecanismos de deformación (STZ, nucleación, coalescencia, fractura)

Escrito el 17/11/2026. Es un análisis en curso: nada pasa a la bitácora hasta que esté maduro. Cada lanzamiento de jobs se propone antes y espera el OK del usuario.

## 0. Preguntas

1. **¿Dónde nacen las STZ?** ¿La estructura local inicial (icosaedros de Cu, Frank-Kasper de Zr, volumen libre, CN) predice qué regiones se activan primero?
2. **¿Cómo crecen y se unen?** ¿Es nucleación de muchos núcleos que coalescen, o propagación desde unos pocos? ¿A qué ε percola la actividad? ¿Forma una banda o una red?
3. **¿Cambia con la historia térmica?** Comparar 1e12 (unificada), 1e11, 2e10 y qgrad_ext. No se usa qgrad (= 2e10).
4. **¿Qué pasa a deformaciones grandes?** ¿Hay localización, cavitación o fractura? ¿Cuándo y dónde?

## 1. Lo que hay hoy y sus limitaciones

**Hecho** (README, 17/11): tracciones a 300 K con 1e8, 1e9 y 1e10 s⁻¹ en x, y, z; a 450 y 540 K con 1e9. Por caso: Voronoi en 41 cuadros, D²min en 20 ventanas de Δε = 0.01 y η en 7 niveles. Primer resultado: los icosaedros frenan el **inicio** de la plasticidad (P(activo | Cu-ico)/P(activo | otro Cu) ≈ 0–0.4 antes del pico) pero no el flujo (≈ 0.8–1.1 para ε ≥ 0.12).

**Aclaración sobre "una loncha":** los números (fracciones, cocientes, Spearman por celdas, AUC) ya se calculan sobre la **caja entera** (65 536 átomos). La loncha de 12 Å se usa solo en los mapas. El problema estadístico real es otro: **cada número sale de UNA realización**.

| debilidad | por qué importa |
|---|---|
| un vidrio por tasa | la diferencia entre tasas incluye la variación entre muestras, que no conocemos |
| una tracción por dirección | las barras "x, y, z" mezclan anisotropía y ruido térmico de una sola trayectoria |
| ruido térmico | a 300 K, una parte de la actividad en una ventana es azar dinámico, no estructura |
| correlación espacial | átomos vecinos no son independientes: los errores tipo √N subestiman la incertidumbre |
| sin intervalos de confianza | no se puede decir si 0,30 contra 0,40 es una diferencia real |
| cuadros cada Δε = 0,005 | la nucleación y la coalescencia ocurren más rápido que eso |
| caja cúbica periódica de 105 Å | traba la orientación de una banda de corte y no tiene superficies libres para romperse |

## 2. Estadística: robustecer las conclusiones (prioridad 1)

### 2a. Con lo que ya existe (sin jobs)
- **Usar todas las realizaciones disponibles:** x, y, z a 300 K (×3 tasas de deformación) y 450/540 K. Hoy la evolución con ε se hizo solo con x a 1e9.
- **Intervalos de confianza por bootstrap espacial en bloques:** dividir la caja en 4³ = 64 bloques de ~26 Å (más que la longitud de correlación de las STZ, ~10–20 Å), remuestrear bloques y recalcular cada cociente. Así se obtiene un IC 95 % honesto dentro de cada realización.
- **Prueba de permutación** para cada cociente: repartir las etiquetas al azar dentro de la especie, como en `clusters_4tasas.py`. Da el p-valor contra "sin efecto".
- **Combinar realizaciones:** media ± error entre direcciones y tasas de deformación; reportar el tamaño del efecto con su IC, no solo el p-valor.
- **Predicción con validación espacial:** el AUC de `nucleacion.py` ya entrena en una mitad y prueba en la otra. Extenderlo a una validación cruzada por bloques y a otras realizaciones (entrenar en x, probar en y).
- **Entregable:** `datos/estadistica_*.csv` y figuras con barras = IC 95 %. Se responde en qué rango de ε el efecto de los icosaedros es significativo y si difiere entre vidrios.

### 2b. Ensamble isoconfiguracional (jobs; el más informativo)
- La misma estructura inicial, N = 10–20 tracciones con **velocidades iniciales distintas** (Widmer-Cooper y Harrowell). Para cada átomo se define la **propensión** = fracción de réplicas en que es activo antes del pico.
- Separa lo que decide la estructura de lo que es azar térmico. Una correlación entre estructura y propensión es mucho más limpia que contra una sola trayectoria.
- **Diseño:** 4 vidrios × 10 semillas, tracción en x a 300 K y 1e9 solo hasta ε = 0.12 (cubre la nucleación y la primera coalescencia). Dumps cada Δε = 0.0025 hasta 0.12.
- **Costo:** ~170 000 pasos por réplica ≈ 10–12 min → 40 jobs ≈ 8 h-nodo (`standby`). Se puede empezar por 2e10 y 1e12 (20 jobs).
- **Análisis:** mapa de propensión; propensión contra descriptores (ico, FK, volumen, CN, q6, entorno a 6 Å), con R² o AUC; ¿la propensión de los icosaedros cambia con la tasa de enfriamiento?

### 2c. Vidrios independientes (jobs; caro; decidir después de 2a/2b)
- Para afirmar que algo cambia con la tasa de enfriamiento hacen falta ≥ 3 vidrios por tasa, enfriados desde configuraciones de líquido distintas (nueva semilla a 2000 K, con el mismo protocolo).
- **Costo por vidrio:** 1e12 ≈ 11 h (rampa + relajación de 10 ns); 1e11 ≈ 25 h; 2e10 ≈ 90 h. Razonable: +2 vidrios de 1e12 y +2 de 1e11 (~70 h-nodo); 2e10 solo si hace falta.
- Sirve también al resto de la tesis (barras de error entre muestras en Tg, estructura, κ).

## 3. Nucleación, crecimiento y coalescencia de STZ (prioridad 2)

### 3a. Ventana fina alrededor del pico (jobs cortos)
- Desde los restarts de ε = 0.05 (y desde ε = 0 si hace falta) rehacer hasta ε = 0.12 con dumps cada **Δε = 0.001** (1 ps a 1e9 s⁻¹), en x, y, z de los 4 vidrios: 12 jobs de ~5–10 min. Si se hace 2b, puede salir de ahí.
- **Seguimiento de clusters entre ventanas:** cada cluster activo (D²min > umbral) se une con los de la ventana siguiente si comparten átomos (solapamiento ≥ 1 átomo o Jaccard > umbral). Eventos: **nacimiento**, **crecimiento**, **fusión** (coalescencia), **desaparición**. Medir:
  - tasa de nacimiento contra ε;
  - distribución de tamaños (¿ley de potencia, avalanchas?);
  - fracción de crecimiento por fusión contra nucleación nueva;
  - ε de percolación (cluster que cruza la caja);
  - ¿los núcleos nuevos aparecen lejos o cerca de los existentes?
- **Relación con la estructura:** propiedades del entorno (ico/FK en r < 6–10 Å) de los sitios donde nace un cluster, contra los sitios al azar y contra los sitios donde un cluster crece.

### 3b. Seguir deformando hasta grandes ε y fractura (jobs)
- **Paso 1 (barato):** continuar las tracciones de 1e9 en x de los 4 vidrios desde el restart de ε = 0.20 hasta ε = 1.0 o hasta que σ → 0. Mismo `fix deform` y NPT lateral; input nuevo `in.tension_cont` que lee el restart y fija el mismo L0 del log original. Dumps cada Δε = 0.01. ~800 000 pasos ≈ 1 h por caso.
- **Qué medir:**
  - curva σ–ε completa;
  - densidad local y **poros** (OVITO `ConstructSurfaceModifier` o celdas de Voronoi grandes): ε y lugar de la primera cavidad, crecimiento y coalescencia de poros;
  - localización (perfil de η a lo largo del eje de tracción, cuello);
  - T local en las zonas activas.
- **Cuidado:** con condiciones periódicas no hay superficies libres. "Romper" = un poro que atraviesa toda la sección y σ → 0. La sección se afina como ~1/(1+ε).
- **Paso 2 (si el paso 1 lo justifica):** geometrías que permitan una banda de corte:
  - caja alargada (aspecto 3:1 en el eje de tracción, `replicate` → ~200 000 átomos) para que entre una banda a ~45°;
  - probeta con superficies libres en y, z (tipo nanopilar; requiere relajar las superficies).

  Costo ~3–5 h por ensayo; se diseña con el OK del usuario.
- **Prueba complementaria (opcional):** deformación uniaxial *sin* contracción lateral (alta triaxialidad) para forzar cavitación y medir la resistencia a la cavitación.

## 4. Visualización 3D (en paralelo, sin jobs)

Probado el 17/11: OVITO 3.15 renderiza sin pantalla con **TachyonRenderer** (caja entera, 400×400, 3.4 s). En anaconda hay **plotly 5.9** y **scikit-image 0.22**.

1. **Renders 3D estáticos (OVITO + Tachyon):** mostrar solo lo relevante (átomos activos coloreados por cluster o por ventana de nacimiento, clusters de icosaedros semitransparentes, el resto oculto o como contorno de la caja). Varias vistas.
2. **Animaciones contra ε:** un cuadro por ventana (nucleación → coalescencia → percolación), armadas como GIF o MP4.
3. **Visor interactivo HTML (plotly):** puntos 3D de los átomos activos (2–6 mil, liviano) y **isosuperficies** de los campos suavizados (η o D²min promediados en una grilla de ~3 Å, `marching_cubes` de scikit-image → `Mesh3d`), más los clusters de icosaedros. Un control deslizante para ε. Se abre en el navegador sin instalar nada; también se puede publicar como página privada.
4. **Para OVITO de escritorio:** exportar dumps livianos (id, tipo, posición, η, D²min, ico, cluster) para abrirlos en la computadora del usuario.

Los mapas 2D quedan solo como apoyo; los cocientes y las curvas contra ε llevan el mensaje cuantitativo.

## 5. Orden propuesto y puntos de decisión

| etapa | qué | jobs | costo | decide el usuario |
|---|---|---|---|---|
| A | 2a (estadística con lo existente) + 4.1/4.3 (3D con lo existente) | no | horas de análisis | qué figuras 3D le sirven |
| B | 3a ventana fina (12 jobs) + 3b paso 1 (4 jobs) | sí | ~6 h-nodo | OK antes de lanzar |
| C | 2b ensamble isoconfiguracional (20–40 jobs) | sí | 4–8 h-nodo | OK, y cuántas semillas |
| D | 3b paso 2 (caja alargada o nanopilar) | sí | ~20–40 h-nodo | solo si B lo justifica |
| E | 2c vidrios independientes | sí | ~70 h-nodo | solo si C deja dudas sobre el efecto de la tasa |

Hace falta que corra la tracción de la 1e12u (1e9, x) para tener los 4 vidrios en A.

## 6. Criterios para considerar "maduro" el análisis (y recién ahí, la bitácora)

- Cada afirmación con su IC 95 % (bootstrap por bloques) y la prueba contra el azar.
- Los resultados se repiten en x, y, z, y en el ensamble isoconfiguracional.
- Las diferencias entre vidrios superan la variación entre realizaciones (y entre vidrios independientes, si se hace E).
- Una figura 3D que muestre nucleación → coalescencia sin ambigüedad.
