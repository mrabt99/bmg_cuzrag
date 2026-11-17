# transporte/: viscosidad y conductividad térmica (Green-Kubo y NEMD)

Viscosidad η y conductividad térmica κ del vidrio y del líquido Zr47.5Cu47.5Ag5 (potencial EAM) con Green-Kubo, y κ también por NEMD. Cada subdirectorio con README explica su método y resultados; los números finales están en los `.txt` y `.csv` que escribe cada script.

## Qué hay en el directorio

| ruta | qué es |
|---|---|
| `cal_<T>K_<tasa>/` | una temperatura y una tasa: NPT 200 ps + 3 réplicas NVE de 1 ns (presión cada 2 fs, J cada 5 fs, MSD). Se crean con `make_cal_dir.sh` |
| `cal_800K_2e10/`, `cal_900K_2e10/` | además, corridas largas para η (`in.nve_long`, tramos de 1 ns) |
| `size_8192/` | vidrio de 8192 átomos (como el paper) y prueba de tamaño de κ; ver su README |
| `nemd/` | κ por NEMD (`fix ehex`), 3 largos de celda y prueba de linealidad; ver su README |
| `vdos/` | VDOS por especie del vidrio 2e10 a 300 K (VACF, 20 orígenes de 4 ps; `in.vdos`, `analiza_vdos.py` → `vdos_300K_2e10.csv`). Figuras en `bitacora/figuras/fig_tr_vdos.py` y `fig_tr_kappa_frecuencia.py` (Re κ(ν)) |
| `verif_argon/`, `verif_acf/` | validación con el ejemplo de argón de LAMMPS y de `fix ave/correlate` |
| `eta_T.py`, `eta_vft.py`, `eta_largo.py` | η(T), κ(T) y D(T) en el líquido; ajuste VFT; corridas largas de 800 y 900 K |
| `kappa_vidrio.py` | κ en el vidrio con el criterio de ventana (trapecio, 0.5-2 ps, test de pendiente) |
| `verif_integridad.py` | chequea que las salidas `cal_r*_*.dat` no estén corruptas (p. ej. tras un relanzamiento) |
| `lanzar_fase2_2e10.sh` | lanzador de la Fase 2 (con guarda contra relanzamientos) |

Reglas prácticas: integrar con trapecio; no restar la media en la ACF; jobs de menos de 4 h con `-q standby`; `srun --mpi=pmi2`; toda corrida a alta T con `fix halt` por distancia mínima (< 1.5 Å) y densidad (< 6.0 g/cm³).

## Posibles extensiones de la tesis (anotado 2026-11-07)

### 1. Contribución de cada tipo de modo vibracional: difusones, propagones y locones

**Motivación.**
- Pregunta abierta: el aporte de locones y propagones a κ.
- En la clasificación de Allen y Feldman (Phys. Rev. B 48, 12581, 1993; Philos. Mag. B 79, 1715, 1999), los modos de un amorfo son de tres tipos:
  - **propagones:** baja frecuencia, casi ondas planas, camino libre largo;
  - **difusones:** el grueso del espectro, extendidos pero no ondulatorios; transportan calor por acoplamiento entre modos;
  - **locones:** alta frecuencia, localizados; casi no transportan.

**Lo que ya sabemos sin descomponer.**
- κ por NEMD casi no cambia entre L = 210 y 838 Å (0.695 → 0.722, no significativo con errores realistas), y GK no cambia entre 8192 y 65536 átomos.
- El aporte de los propagones que cabe en estas cajas es, entonces, chico (≲10 %); con estos datos no se puede medir.
- La cola de la ACF entre 2 y 20 ps es compatible con cero.

**Cómo hacerlo** (el LAMMPS de Negishi tiene el paquete PHONON: `dynamical_matrix`, `third_order`, `fix phonon`):
1. **Densidad de estados vibracionales (VDOS).** Transformada de Fourier de la autocorrelación de velocidades, por especie. Barato: unos ps de NVE a 300 K guardando velocidades. **Hecho 2026-11-14 (2e10, 300 K, `vdos/`, job 45840090, 5 min):** banda hasta ~8 THz; máximos Ag 2,3, Zr 3,3, Cu 3,5 THz (Cu con hombro a ~5,7 THz); total 3,4 THz, ⟨ν⟩ = 4,0 THz. ν_D efectiva de g ≈ 3ν²/ν_D³ a 0,4-1,2 THz: 3,4-3,6 THz, contra 4,4 THz de las elásticas (exceso de modos blandos, aunque la resolución de la ventana de 3 ps también suma). **Revisado el mismo día:** con la ventana entre 1,5 y 3,9 ps, g/ν² dividido por el nivel de Debye (3/4,36³ = 0,0362 THz⁻³) da 2,2 / 2,05 / 1,85 / 1,3 / 1,0 a 0,5 / 0,7 / 0,9 / 1,5 / ~2 THz, igual para todas las ventanas. El exceso de ~2× a 0,5-1 THz es real, no un efecto de la resolución, que solo afecta por debajo de ~0,4 THz. Pero g/ν² crece en forma monótona hacia ν → 0 hasta el límite de resolución: no hay máximo resuelto, así que **no se puede afirmar que sea un pico bosónico**. Para resolverlo: VACF con tramos más largos (resolución < 0,1 THz); el modo transversal más bajo de la caja está en v_T/L ≈ 0,15 THz. **Hecho 2026-11-15** (`run_vdos_largo.sbatch`, 20 tramos de 20 ps, job 45842063, 23 min; `analiza_reducida.py v300_2e10_L20 3 5 10 15 19` → `reducida_v300_2e10_L20.csv`; figura `bitacora/figs/tr_vdos_reducida.png`). Con t_c = 15 ps (resolución 0,05 THz), g/ν² es una meseta de 2,0-2,2× Debye entre 0,15 y 0,7 THz y baja por encima (1,75 a 1 THz, 1,0 a 2 THz), igual con t_c = 10-19 ps. La subida hacia ν → 0 de la corrida de 4 ps era la ventana. Los máximos locales (0,15, 0,22, 0,35 THz) coinciden con los modos transversales de la caja (|n|² = 1, 2, 5), lo que de paso confirma v_T = 1611 m/s y, con eso, el nivel de Debye. Las ondas planas de la caja por debajo de 0,6 THz son 518, contra 512 de Debye: la discretización no explica el factor 2 (~500 modos de más que no son ondas planas). **Exceso real, sin máximo resuelto:** si hay pico bosónico, está por debajo de ~0,5 THz y es muy ancho. Para bajar de 0,15 THz hace falta una caja 2×2×2 (L = 210 Å, primer modo 0,077 THz, ~8× el costo). **Caja 2×2×2 hecha 2026-11-16** (`run_vdos_rep2.sbatch`: `in.vdos` con `-var rep 2`, mismo vidrio replicado, 524288 átomos, velocidades nuevas con semilla 6006, 10 tramos de 40 ps; job 45925782, 73 min en 128 núcleos; T = 300,2 K, E_total sin deriva; `analiza_reducida.py v300_2e10_rep2 5 10 20 30 38`). Con t_c = 30 ps (resolución 0,024 THz), g/ν² ÷ Debye = 2,20 / 2,03 / 2,02 / 2,09 en 0,06-0,1 / 0,1-0,15 / 0,15-0,3 / 0,3-0,6 THz. Coincide con la caja chica donde se solapan, así que no hay efecto de tamaño. Los máximos locales vuelven a caer en los modos transversales de la caja (0,077, 0,109, 0,172-0,188 THz), y las ondas planas bajo 0,3 THz son 518 contra 512 de Debye. **El exceso de ~2× es plano desde 0,07 hasta 0,7 THz, sin máximo: no hay pico bosónico resoluble.** Como en ν → 0 g/ν² tiene que volver al nivel de Debye, falta entender por qué no lo hace ni a 0,07 THz. Prueba propuesta, no hecha: la misma VDOS a 50 K (anarmonicidad casi nula) contra el Debye de las Cij a 0 K (ΘD = 228 K, ν_D = 4,75 THz). Re κ(ν) ocupa la misma banda (también se anula hacia 8 THz) pero está corrida a frecuencias más bajas (máximo 2,5 THz) y es finita en ν→0 (0,82), mientras que g → 0. Ojo: J es cuadrático en las vibraciones, así que κ(ν) refleja pares de modos (frecuencias de batido) y no "el aporte de los modos de frecuencia ν". Para separar por modo hacen falta los pasos 2-4.
2. **Modos armónicos.** Minimizar el vidrio de 8192 átomos (estructura inherente), calcular la matriz dinámica con `dynamical_matrix` (3N = 24 576) y diagonalizarla.
   - La matriz ocupa ~4.8 GB en doble precisión; la diagonalización es factible en un nodo.
   - Por modo se obtienen la **razón de participación** (identifica a los locones) y el **carácter ondulatorio** (factor de estructura dinámico del autovector), que separa propagones de difusones (criterio de Ioffe-Regel).
3. **κ armónica de Allen-Feldman** (difusividad de cada modo a partir de la matriz dinámica y el operador de flujo de calor). Comparada con GK/NEMD separa la parte armónica de la anarmónica.
4. **Descomposición espectral del flujo en NEMD** (flujo de calor espectral, Sääskilahti et al., Phys. Rev. B 90, 134312, 2014). Da κ(ω) directamente de las corridas de `nemd/`, guardando velocidades y fuerzas entre bloques vecinos.

**Resultado esperado.** Fracción de κ por tipo de modo, y una respuesta cuantitativa a la pregunta de los propagones.

### 2. Contribución electrónica a κ

**Motivación.**
- El potencial EAM no tiene electrones: tanto nuestro κ como el del paper son solo la parte vibracional (iónica).
- En un vidrio metálico real buena parte del calor la llevan los electrones. Para comparar con un κ experimental hay que sumarla.

**Opciones, de más simple a más costosa:**
1. **Ley de Wiedemann-Franz:** κ_e = L₀ σ T, con L₀ = 2.44 × 10⁻⁸ W·Ω/K². La conductividad eléctrica σ(T) se toma de mediciones de resistividad de vidrios Zr-Cu(-Ag) en la literatura (hay que buscarlas). Es inmediato y da κ_total = κ_e + κ_vib para comparar con experimentos. Limitación: en sistemas muy desordenados L puede apartarse de L₀; conviene citar la incertidumbre.
2. **Primeros principios (DFT + Kubo-Greenwood)** sobre configuraciones de MD de pocos cientos de átomos: σ y κ_e sin datos experimentales. Requiere un código DFT (VASP o Quantum ESPRESSO) y bastante cómputo.
3. **Acople electrón-fonón** (modelo de dos temperaturas): más relevante para irradiación que para κ en equilibrio. Probablemente queda fuera.

**Agregado 2026-11-14 (para hacer en el corto plazo).**
- **Alcance.** Solo κ (y σ, Seebeck) necesita electrones. η y D son coeficientes iónicos: los electrones entran solo vía el potencial (adiabático), que el EAM ya incluye, y su aporte directo a η es despreciable. Lo de GK/NEMD de este directorio es κ de red (κ_vib).
- **Orden de magnitud (Wiedemann-Franz).** Vidrios Cu-Zr: ρ ≈ 150-200 μΩ·cm (a confirmar para Zr-Cu-Ag). Con ρ = 180 μΩ·cm: κ_e(300 K) ≈ 2.44e-8 × 300 / 1.8e-6 ≈ 4 W/m/K, unas 5 veces κ_vib ≈ 0.8. Los electrones llevan ~80-85 % del calor. Por la regla de Mooij (ρ alta ⇒ ρ casi independiente de T, coeficiente de temperatura chico y negativo), κ_e crece aproximadamente lineal con T, mientras que κ_vib es plana en el vidrio.
- **Opción nueva: teoría de Ziman / Faber-Ziman extendida.** Resistividad de metales amorfos y líquidos a partir de los factores de estructura parciales S_ij(q) más los desfasajes de dispersión electrón-ion (t-matriz de muffin-tin). Las S_ij(q) salen de las g_ij(r) parciales de la Etapa 4 (`structural_analysis/etapa4/`). La parte electrónica es un cálculo atómico chico, sin DFT de la celda. Es semicuantitativa (±20-30 % en ρ), pero da la dependencia con T y la regla de Mooij y conecta σ con la estructura y con la tasa de enfriamiento.
- **Sobre la opción 2 (DFT + Kubo-Greenwood).** Celdas de ~150-250 átomos cortadas de nuestros vidrios (relajación DFT corta), varias configuraciones. Da σ(ω) y los coeficientes de Onsager L₁₁, L₁₂, L₂₂ → κ_e sin suponer Wiedemann-Franz, el número de Lorenz real y el Seebeck. Caro por los estados 4d del Zr y los puntos k, pero factible en Negishi.
- **Corrección a la opción 3.** El modelo de dos temperaturas (`fix ttm`) recibe κ_e como dato de entrada; no lo calcula. No sirve para obtener κ_e.
- **Plan propuesto.**
  1. Búsqueda bibliográfica de ρ(T) y κ medidos en vidrios Zr-Cu-Ag / Cu-Zr.
  2. Wiedemann-Franz → tabla κ_e(T) y κ_total = κ_e + κ_vib (300-600 K), con la incertidumbre de L₀.
  3. En la tesis, decir explícitamente que la MD da solo κ_vib. La comparación con un experimento (κ ≈ 4-6 W/m/K, casi todo electrónico) pesa más que la discrepancia 0.8 contra 2.4 con Méndez.
  4. Opcional, como aporte propio: Ziman con nuestras S_ij(q), comparando tasas. Si hay tiempo y cómputo, Kubo-Greenwood en 2-3 configuraciones como validación.
