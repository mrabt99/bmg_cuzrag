# Prueba de la hipótesis del término convectivo

**Pregunta:** κ es plana en el vidrio (~0.79 W/m/K) y sube en el líquido (+0.021 por 100 K, hasta 0.88 a 1200 K). ¿La subida viene del término convectivo del flujo de calor, que crece con la difusión, y no de la conducción por vibraciones?

## Flujo de calor y su partición

`compute heat/flux` da J = Σ eᵢvᵢ − Σ Sᵢvᵢ: la parte convectiva (componentes 4-6) y la de virial (la diferencia).

En una mezcla, la parte convectiva depende del cero de energía de cada especie. Con momento total nulo, Σ ē vᵢ = 0 para una sola especie, pero con tres, Σₛ ēₛ Jₛ (con Jₛ = Nₛ v_cm,s) no se anula: la interdifusión transporta la energía media de cada especie. Lo usual es definir el flujo de calor restando las entalpías parciales, J' = J − Σₛ hₛ Jₛ. Esa parte crece con la difusión y es candidata a explicar la subida de κ en el líquido.

## Qué se guarda (`in.nve_conv`)

Es igual a `cal_1100K_2e10/in.nve`; solo agrega salidas, cada 5 fs como J:
- `conv_r<N>_conv.dat`: convectiva total (hf[4-6]) y por especie (hZr, hCu, hAg [4-6]), en eV·Å/ps.
- `conv_r<N>_vcm.dat`: velocidad del centro de masa de Zr, Cu y Ag (Å/ps).
- `conv_r<N>_hflux.dat`, `_press.dat`, `_msd.dat` y `_thermo.dat`, con el mismo formato que `cal_*`. Así, estas réplicas además se suman a las existentes.

## Corridas

| T (K) | partida | réplicas | semillas |
|---|---|---|---|
| 600 | `../cal_600K_2e10/cal_npt.restart` (vidrio, control) | 3 × 1 ns | 7007, 8008, 9009 |
| 900 | `../cal_900K_2e10/cal_npt.restart` | 3 × 1 ns | 7007, 8008, 9009 |
| 1200 | `../cal_1200K_2e10/cal_npt.restart` | 3 × 1 ns | 7007, 8008, 9009 |

NVT 20 ps + NVE 1 ns, 48 tareas, standby, ~1 h por job: 9 jobs, sin NPT (se reusa el de `cal_*`, solo lectura).

`./lanzar_conv.sh` encola las 9 corridas (o `./lanzar_conv.sh 900` para una sola T). Saltea las réplicas que ya tengan salida, log o job en cola. **Lanzar solo con go-ahead.**

## Prueba de sintaxis (2026-11-11)

900 K, 8 procesos en el login, NVT de 100 pasos + NVE de 200, en el scratchpad: rc = 0, sin errores.
- La suma de la convectiva por especie es igual a la total: en la primera fila, en x, −972.5 + 2093.7 − 507.2 = 614.0 contra 614.0; lo mismo en y y z.
- Σ mₛNₛv_cm,s ≈ 0: momento total nulo, dentro del redondeo.
- Cada parte por especie está dominada por ēₛNₛv_cm,s: para el Zr en x, ≈ −910 contra −972 medido. Es una sola muestra, sin valor estadístico.

## Resultado (2026-11-11, jobs 45550794-802, `analiza_conv.py` → `analiza_conv.txt/png`)

Las 9 corridas terminaron COMPLETED: sin errores ni halts, 200 001 filas alineadas, suma por especie igual al total (10⁻⁴ relativo, redondeo de la salida), momento total ≈ 0 (8×10⁻⁴ relativo, redondeo), E_total conservada, T dentro de ±3 K.

κ en la ventana 0.5-2 ps (W/m/K):

| T (K) | total | virial-virial | conv-conv | cruzada | **proyectado** (flujo de masa nulo) |
|---|---|---|---|---|---|
| 600 | 0.771 ± 0.017 | 0.687 ± 0.015 | 0.010 | +0.074 | = total (sin difusión; proyección mal condicionada) |
| 900 | 0.822 ± 0.022 | 0.686 ± 0.017 | 0.053 | +0.082 | **0.788 ± 0.021** |
| 1200 | 0.925 ± 0.018 | 0.696 ± 0.013 | 0.110 | +0.118 | **0.822 ± 0.014** |

- **La parte de vibraciones (virial-virial) es plana**: 0.687 → 0.686 → 0.696. De la subida del total entre 600 y 1200 K (+0.154), 0.100 es convectiva-convectiva, 0.044 cruzada y solo 0.009 virial. **Hipótesis confirmada: la subida de κ en el líquido es convectiva.**
- **El κ total del líquido depende del cero de energía**: con el cero corrido c = −5…+5 eV, 0.818-0.883 a 900 K y 0.919-1.061 a 1200 K; a 600 K, 0.771-0.774. Con tres especies, Σᵢvᵢ ≠ 0 aunque el momento total sea nulo, y la interdifusión transporta la energía de cada especie. El κ total de Green-Kubo de la mezcla en el líquido no está bien definido; el bien definido es el proyectado, a flujo de masa nulo (Onsager, L_qq − L_qd L_dd⁻¹ L_dq), que no depende de c y es el que mide un experimento o NEMD estacionario.
- Energía efectiva por especie (pendiente de la convectiva de la especie contra su flujo): Zr −6.0, Cu −3.5, Ag −2.7 eV, parecida a las energías de cohesión de los elementos puros (≈ 6.3, 3.5 y 2.9 eV). Es un control de coherencia.
- **κ proyectado**: sube solo de ~0.77-0.80 (vidrio) a 0.79 (900 K) y 0.82 (1200 K). Es estable entre ventanas: 900 K 0.788 / 0.794 / 0.828 / 0.811 (0.5-2, 2-5, 5-10, 10-20 ps); 1200 K 0.822 / 0.819 / 0.777 / 0.710. La caída de 1200 K después de 5 ps acompaña a la del virial-virial (ruido de la cola). Las integrales de interdifusión convergen a 1-2 ps a 1200 K; a 900 K tienen una cola negativa (−20 % entre 2 y 10 ps) que no cambia el κ proyectado dentro del error.
- El pico de κ(t) a ~0.1 ps en el líquido es casi todo convectivo; el κ proyectado no tiene pico.
- 600 K con semillas nuevas: total 0.771 ± 0.017, contra 0.825 ± 0.019 de r1-r3 (−0.054, 2.1σ). La dispersión entre juegos de réplicas es mayor que el sem de bloques en ese punto (que además no pasa el test de meseta).

**Implicancias:** ver el plan (`plan_transporte_GK.md`, sección "Hipótesis del término convectivo").

## Análisis previsto (escrito antes de tener los datos; se hizo en `analiza_conv.py`)

1. κ separado en convectiva-convectiva, virial-virial y cruzada, con el mismo criterio que `kappa_tasas.py` (trapecio, bloques de 100 ps, ventana 0.5-2 ps, T medida).
2. κ con J' = J − Σₛ ēₛ Jₛ, usando la energía media por átomo de cada especie (de la convectiva y v_cm). Si la subida entre 600 y 1200 K desaparece con J', la subida es interdifusión y no conducción por vibraciones.
3. Sensibilidad al cero de energía: κ con eᵢ → eᵢ − c para distintos c. El κ total no debería depender de c; si depende, la partición no está bien definida.

Cuidado: la energía media por átomo de cada especie no es exactamente la entalpía parcial. A P = 0 la diferencia es la dependencia con la composición, que acá no se puede medir. Si el resultado depende de esa elección, hay que reportarlo como límite.
