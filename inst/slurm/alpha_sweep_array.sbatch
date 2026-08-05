#!/bin/bash
#SBATCH --job-name=smfsim_alpha_grid
#SBATCH --account=EDIT_ME_ACCOUNT          # <-- EDIT ME
#SBATCH --qos=EDIT_ME_QOS                  # <-- EDIT ME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jobrant@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48gb
#SBATCH --time=08:00:00
#SBATCH --array=1-15
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# NOTE: these paths are relative to the directory you run `sbatch` FROM, not to
# $WORKDIR, and SLURM opens them BEFORE this script runs -- so the `mkdir -p
# logs` below cannot help. Run `mkdir -p logs` in your submit directory first,
# and submit from the same place every time, or the output silently goes
# somewhere you are not looking (or is lost).

# One GRID CELL per array task (a within_alpha x between_alpha pair).
#
# --array=1-15 matches the default 5 x 3 grid in run_alpha_sweep.R. If you
# change WITHIN_GRID or BETWEEN_GRID, change this range to match; an
# out-of-range task id fails fast with an explicit error.
#
# Add %5 (i.e. --array=1-15%5) to cap concurrency if the cluster is busy or if
# 15 tasks reading the same allc files causes I/O contention.
#
# The script detects SLURM_ARRAY_TASK_ID and runs ONLY that cell, writing to
# its own w<within>_b<between>/ directory. It deliberately does NOT write the
# combined CSVs in array mode, to avoid write races -- that is the stitch job.
#
# Time: 8 h suits FAST_TUNE <- TRUE (chr1). Genome-wide needs far more; confirm
# only the single chosen pair genome-wide rather than the whole grid.

set -euo pipefail

WORKDIR=EDIT_ME_PATH_TO_SMFsim             # <-- EDIT ME (package root)
module load EDIT_ME_R_MODULE               # <-- EDIT ME, e.g. R/4.4

cd "$WORKDIR"
mkdir -p logs

echo "host        : $(hostname)"
echo "array job   : ${SLURM_ARRAY_JOB_ID}"
echo "array task  : ${SLURM_ARRAY_TASK_ID}"
echo "workdir     : $(pwd)"
echo "started     : $(date)"

Rscript inst/scripts/run_alpha_sweep.R

echo "finished    : $(date)"
