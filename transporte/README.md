# transporte/: viscosidad y conductividad térmica (Green-Kubo y NEMD)

Viscosidad η y conductividad térmica κ del vidrio y del líquido Zr47.5Cu47.5Ag5 (potencial EAM) con Green-Kubo, y κ también por NEMD. El plan, los resultados y todas las decisiones de método están en **`plan_transporte_GK.md`**, que es la fuente de verdad, con bloques fechados. Este README solo orienta.

## Qué hay en el directorio

| ruta | qué es |
|---|---|
| `cal_<T>K_<tasa>/` | una temperatura y una tasa: NPT 200 ps + 3 réplicas NVE de 1 ns (presión cada 2 fs, J cada 5 fs, MSD). Se crean con `make_cal_dir.sh` |
| `cal_800K_2e10/`, `cal_900K_2e10/` | además, corridas largas para η (`in.nve_long`, tramos de 1 ns) |
| `size_8192/` | vidrio de 8192 átomos (como el paper) y prueba de tamaño de κ; ver su README |
| `nemd/` | κ por NEMD (`fix ehex`), 3 largos de celda y prueba de linealidad; ver su README |
| `verif_argon/`, `verif_acf/` | validación con el ejemplo de argón de LAMMPS y de `fix ave/correlate` |
| `mail_seba/`, `mail_seba2/` | mails a Jaroszewicz y Méndez, figuras, tarballs y el notebook y los datos de Seba |
| `eta_T.py`, `eta_vft.py`, `eta_largo.py` | η(T), κ(T) y D(T) en el líquido; ajuste VFT; corridas largas de 800 y 900 K |
| `kappa_vidrio.py`, `mail_seba/kappa_vs_T.py` | κ en el vidrio con el criterio de ventana (trapecio, 0.5-2 ps, test de pendiente) |
| `replica_metodo_seba.py` | aplica el procedimiento del notebook de Seba a nuestras series |
| `verif_integridad.py` | chequea que las salidas `cal_r*_*.dat` no estén corruptas (p. ej. tras un relanzamiento) |
| `lanzar_fase2_2e10.sh` | lanzador de la Fase 2 (con guarda contra relanzamientos) |

Reglas prácticas: integrar con trapecio; no restar la media en la ACF; jobs de menos de 4 h con `-q standby`; `srun --mpi=pmi2`; toda corrida a alta T con `fix halt` por distancia mínima (< 1.5 Å) y densidad (< 6.0 g/cm³).

## Posibles extensiones de la tesis (anotado 2026-11-07)

### 1. Contribución de cada tipo de modo vibracional: difusones, propagones y locones

**Motivación.**
- Jaroszewicz preguntó por el aporte de locones y propagones a κ.
- En la clasificación de Allen y Feldman (Phys. Rev. B 48, 12581, 1993; Philos. Mag. B 79, 1715, 1999), los modos de un amorfo son de tres tipos:
  - **propagones:** baja frecuencia, casi ondas planas, camino libre largo;
  - **difusones:** el grueso del espectro, extendidos pero no ondulatorios; transportan calor por acoplamiento entre modos;
  - **locones:** alta frecuencia, localizados; casi no transportan.

**Lo que ya sabemos sin descomponer.**
- κ por NEMD casi no cambia entre L = 210 y 838 Å (0.695 → 0.722, no significativo con errores realistas), y GK no cambia entre 8192 y 65536 átomos.
- El aporte de los propagones que cabe en estas cajas es, entonces, chico (≲10 %); con estos datos no se puede medir.
- La cola de la ACF entre 2 y 20 ps es compatible con cero.

**Cómo hacerlo** (el LAMMPS de Negishi tiene el paquete PHONON: `dynamical_matrix`, `third_order`, `fix phonon`):
1. **Densidad de estados vibracionales (VDOS).** Transformada de Fourier de la autocorrelación de velocidades, por especie. Barato: unos ps de NVE a 300 K guardando velocidades.
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
