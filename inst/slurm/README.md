# SLURM templates

Submit from the **package root** (the directory containing `DESCRIPTION`), not
from this directory — every script `cd`s to `$WORKDIR` and the R scripts call
`devtools::load_all(".")`.

Fill in the four `EDIT ME` values at the top of each file: account, qos,
workdir, and the R module name.

## Order

```bash
mkdir -p logs

# 1. Calibration gate. Do NOT skip - if the null is dirty, the sweep is waste.
sbatch inst/slurm/validate_null.sbatch

# 2. The grid: one cell per array task.
ARRAY_JOB=$(sbatch --parsable inst/slurm/alpha_sweep_array.sbatch)

# 3. Stitch + plots + separability report, after every cell succeeds.
sbatch --dependency=afterok:$ARRAY_JOB inst/slurm/alpha_sweep_stitch.sbatch
```

## Array size must match the grid

`--array=1-15` assumes the default 5 x 3 grid in `run_alpha_sweep.R`
(`WITHIN_GRID` x `BETWEEN_GRID`). If you change either grid, change the array
range to match. To check the cell count without running anything:

```bash
Rscript -e 'length(seq(0.1,0.9,0.2)) * length(c(0.5,0.7,0.9))'
```

An out-of-range task id fails fast with an explicit error rather than running
the wrong cell.

## Resume

Both scripts resume per cell: a cell whose `null_simulation_results.csv` and
`spikein_simulation_results.csv` already exist is read from cache. So a
partially-failed array can simply be resubmitted — finished cells are skipped.
To force a cell to re-run, delete its `w<within>_b<between>/` directory.

## Concurrency

`--array=1-15%5` caps it at 5 simultaneous tasks. Worth using if the cluster is
busy or if all 15 tasks reading the same `allc` files causes I/O contention.
