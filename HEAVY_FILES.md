# Archivos pesados excluidos de git

Este repo excluye (vía `.gitignore`) todo archivo individual mayor a 20MB — en la
práctica, son las trayectorias MD crudas (`.dump`), que suman ~2.9GB. Todo lo demás
(scripts, potencial EAM, restarts/checkpoints, logs, resultados de BOO, plots, CSVs)
sí está en el repo.

Esta lista categoriza esos archivos pesados por prioridad de backup manual (scp/rsync
desde Negishi a tu PC), para el caso de perder acceso al cluster.

## Prioridad alta — no respaldados en ningún otro lado (~1.67GB)

Trayectorias de producción de las corridas propias (1e11 y la referencia 1e12 del
profesor). Sin esto, para recuperarlas habría que re-rodar los segmentos de MD
(horas de cómputo HPC por segmento) o, en el caso de la 1e12, no serían recuperables
en absoluto porque es la corrida original que compartió el profesor.

- **`rate_1e12_reference/quench_N65536_1e12.dump`** (753MB) + `relax300K_N65536_1e12.dump` (41MB)
  — corrida original del profesor a 1e12 K/s. **Irreemplazable.**
- **`rate_1e11/prod_1e11/quench_N65536_1e11.dump.seg1..seg16`** (~48MB c/u, ~815MB) +
  `relax300K_N65536_1e11.dump` (40MB) + `sweep_1e11/quench_N65536_1e11.dump.seg1_global` (48MB)
  — trayectoria completa de producción a 1e11 K/s, 17 segmentos de 100K.

```bash
# Ejemplo para bajar todo este grupo con rsync (correr desde tu PC):
rsync -avz aaabt@negishi.rcac.purdue.edu:/scratch/negishi/aaabt/bmg_cuzrag/rate_1e12_reference/*.dump ./rate_1e12_reference/
rsync -avz aaabt@negishi.rcac.purdue.edu:/scratch/negishi/aaabt/bmg_cuzrag/rate_1e11/prod_1e11/*.dump ./rate_1e11/prod_1e11/
rsync -avz aaabt@negishi.rcac.purdue.edu:/scratch/negishi/aaabt/bmg_cuzrag/rate_1e11/sweep_1e11/*.dump ./rate_1e11/sweep_1e11/
```

## Prioridad media — ya respaldados en GitHub, pero en otro repo (~857MB)

- **`rate_2e10/prod_2e10/quench_N65536_2e10.dump.seg1..seg16`** + `relax300K_N65536_2e10.dump`
  + `sweep_2e10/quench_N65536_2e10.dump.seg1_global`
  — trayectoria completa a 2e10 K/s. Esto **ya está pusheado** en el repo viejo
  `github.com/mrabt99/meltquench2e10` (commit `8026487`), que se dejó intacto al
  armar este repo nuevo. No es urgente bajarlo aparte, pero tampoco está en tu PC
  todavía — si querés todo en un solo lugar, cloná ese repo o traelo con rsync igual
  que arriba.
- El `.git` viejo de ese repo (con toda esta data ya empaquetada, ~725MB) quedó
  además como backup local en `bmg_cuzrag/rate_2e10_old_git_backup/` (no se subió,
  está en `.gitignore`).

## Prioridad baja — corrida abandonada, no hace falta bajarla (~433MB)

- `for_local_download/quench_N65536_1e10.dump.validate_seed2` (199MB)
- `rate_2e10/quench_N65536_1e10.dump.validate` (35MB)
- `rate_2e10/quench_N65536_1e10.dump.validate_seed2.part_aa/ab/ac` (199MB, mismo
  archivo partido en 3)

Son diagnósticos de la corrida original a 1e10 K/s que se abandonó (explotó por el
core repulsivo inválido del potencial EAM — ver `rate_2e10/ESTADO_PROYECTO_meltquench_1e10.md`).
No aportan a la metodología final a 2e10 K/s. Se pueden ignorar tranquilamente.
