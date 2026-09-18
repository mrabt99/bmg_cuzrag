# Prompt para el asistente de VSCode — Quench Cu-Zr(-Ag) a 10^10 K/s en Negishi

## Contexto

Tu profesor compartió `/scratch/negishi/aaabt/65536` con dos inputs de LAMMPS:

- `in.quench1e12`: enfría el sistema de 2000 K a 300 K a 10^12 K/s.
- `in.relax300K_1E12` (verificar el nombre exacto del archivo en la carpeta): equilibra a 300 K durante 10 ns.

Los dumps de `in.quench1e12` guardan posiciones y velocidades cada 10 K (170 temperaturas, ~11.2 millones de líneas en total) y hay un restart cada 100 K para poder arrancar corridas a temperaturas intermedias (por ejemplo, curvas tensión-deformación).

## Objetivo

Repetir el mismo protocolo (mismo sistema, mismo potencial EAM, misma T inicial y final, mismo ensemble) pero con una tasa de enfriamiento de 10^10 K/s en lugar de 10^12 K/s. Es 100 veces más lento.

## Instrucciones

1. Antes de tocar nada, leé `in.quench1e12` completo (encabezado incluido) y anotá: potencial usado, N° de átomos/caja, ensemble (NPT/NVT), timestep, Tdamp/Pdamp si los hay, cómo está definida la rampa de enfriamiento, y cómo están calculados los intervalos de dump y restart en pasos.
2. Copiá los inputs a un directorio nuevo, por ejemplo `/scratch/negishi/aaabt/65536_1e10/`. No edites ni borres nada dentro de la carpeta original.
3. Calculá el nuevo número de pasos de quench: con ΔT = 1700 K (2000→300) a 10^10 K/s, el tiempo físico total pasa de ~1.7 ns a ~170 ns. Con el mismo timestep, son ~100 veces más pasos de MD.
4. Ajustá los intervalos de dump (cada 10 K) y restart (cada 100 K) en número de pasos. Como ahora entran ~100 veces más pasos por cada 10 K de descenso, copiar el intervalo original tal cual va a dejar los dumps/restarts en lugares equivocados o generar un volumen de datos muy distinto al esperado. Recalculá el intervalo en pasos a partir de la nueva tasa.
5. Armá el script de SLURM (mismos módulos/versión de LAMMPS que usó tu profesor, si lo podés confirmar) y estimá el walltime antes de pedir el allocation completo.

## Cinturones de seguridad y cuidados para que no explote la corrida

**Antes de lanzar el job grande**

- Corré un test corto (unos pocos miles de pasos, cola interactiva o debug) con el input nuevo antes de mandar el job completo. Confirmá que arranca sin `ERROR`, `Lost atoms`, `Bond/angle/dihedral extent`, ni NaN/Inf en energía, temperatura o presión, y que el volumen de la caja no se dispara.
- Con ese test medí el rendimiento (ns/día) y extrapolá el tiempo real para completar los ~170 ns de quench en la partición que uses. 100 veces más pasos puede no entrar en un solo walltime: si es así, encadená jobs con `sbatch --dependency=afterany`, reiniciando siempre desde el último restart, en vez de pedir un walltime gigante de una sola vez.
- Revisá cuánto espacio en disco vas a necesitar. Si el intervalo de dump quedó mal ajustado (punto 4) podés terminar generando mucho más output del esperado y llenar la cuota del scratch.

**Durante la corrida**

- Guardá restarts con frecuencia razonable (no solo "cada 100 K" si ese tramo ahora dura 100 veces más tiempo real), para no perder horas de cómputo si el job se cae o se corta por walltime.
- Revisá el thermo output (energía potencial, temperatura, presión) a intervalos cortos al principio de la corrida para confirmar que se comporta como el caso original antes de dejarla correr sin supervisión.

**Cuidado específico de este sistema (no solo "explosión" numérica)**

- A tasas tan bajas como 10^10 K/s, Ryltsev et al. (2016, en el proyecto) reportan cristalización espontánea (cristalitos de Cu2Zr) en Cu-Zr a 1.5×10^9 K/s, un orden de magnitud por debajo de la tuya, y señalan que desde 10^10 K/s hacia abajo la estructura ya depende más fuerte de la tasa. No es un problema numérico, pero puede arruinar la muestra como vidrio sin que la corrida "explote" en el sentido de LAMMPS: vigilá una caída abrupta e inusual de la energía potencial durante el enfriamiento (además de la caída suave esperada en la transición vítrea), señal de nucleación cristalina.
- Si el script usa termostato/barostato de Nose-Hoover, confirmá que Tdamp y Pdamp siguen siendo razonables para el mismo timestep. No dependen de la tasa de enfriamiento, así que en principio no hay que tocarlos, pero confirmalo en el input real en vez de asumirlo.

## Al terminar

Reportá si terminó sin errores, cuánto tardó realmente vs. lo estimado, y el aspecto de energía potencial vs. temperatura (inflección clara de transición vítrea, sin caídas raras). Dejá un README corto en la carpeta nueva con la tasa usada, el timestep, el número de pasos y los intervalos de dump/restart en pasos, ya que la carpeta original no tiene uno.
