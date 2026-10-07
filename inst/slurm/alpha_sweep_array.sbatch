#!/bin/bash
#SBATCH --job-name=smfsim_alpha_grid
#SBATCH --account=jobrant          # <-- EDIT ME
#SBATCH --qos=jobrant-b                  # <-- EDIT ME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shenxinyu@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48gb
#SBATCH --time=08:00:00
#SBATCH --array=1-12
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
# --array=1-12 matches the 3 x 4 grid in run_alpha_sweep.R (within 0.3/0.5/0.7
# x between 0.8/0.9/0.95/1.0). If you change WITHIN_GRID or BETWEEN_GRID,
# change this range AND EXPECTED_CELLS in alpha_sweep_stitch to match. An
# out-of-range task id fails, and one failed task stops an afterok-dependent
# stitch job from ever starting -- which is what the old 1-15 range did.
#
# Add %4 (i.e. --array=1-12%4) to cap concurrency if the cluster is busy or if
# 12 tasks reading the same allc files causes I/O contention.
#
# The script detects SLURM_ARRAY_TASK_ID and runs ONLY that cell, writing to
# its own w<within>_b<between>/ directory. It deliberately does NOT write the
# combined CSVs in array mode, to avoid write races -- that is the stitch job.
#
# Time: 8 h suits FAST_TUNE <- TRUE (chr1). Genome-wide needs far more; confirm
# only the single chosen pair genome-wide rather than the whole grid.

set -euo pipefail

module load R/4.5             

mkdir -p logs

echo "host        : $(hostname)"
echo "array job   : ${SLURM_ARRAY_JOB_ID}"
echo "array task  : ${SLURM_ARRAY_TASK_ID}"
echo "workdir     : $(pwd)"
echo "started     : $(date)"

Rscript SMFsim/inst/scripts/run_alpha_sweep.R

echo "finished    : $(date)"
