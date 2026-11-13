#!/bin/bash
#SBATCH --job-name=smfsim_bias_seeds
#SBATCH --account=jobrant          # <-- EDIT ME
#SBATCH --qos=jobrant-b            # <-- EDIT ME
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=96:00:00
#SBATCH --array=1-5
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# The bias experiment (run_bias_experiment.R), once per seed: one array task
# per entry of SEEDS below. Keep --array in step with the number of seeds.
#
# Log paths are relative to the SUBMIT directory and are opened before this
# script runs: run `mkdir -p logs` there first, and submit from the manuscript
# root (the directory holding SMFsim/ and data/).
#
#   sbatch SMFsim/inst/slurm/bias_experiment_seeds.sbatch
#
# Each task writes to  <SMFSIM_SEEDS_OUTPUT>/seed_<seed>/<block>/ . When all
# tasks have finished, combine them:
#
#   Rscript SMFsim/inst/scripts/summarize_seeds.R
#
# Options go on the sbatch line, e.g. the logit efficiency model into its own
# output directory:
#   sbatch --export=ALL,SMFSIM_SEEDS_OUTPUT=results/bias_experiment_seeds_logit,SMFSIM_EXTRA_ARGS='--efficiency_model logit' \
#       SMFsim/inst/slurm/bias_experiment_seeds.sbatch
#
# Each task is a full genome-wide run of every block, the same cost as
# bias_experiment.sbatch.

set -euo pipefail

# === EDIT ME =============================================================
R_MODULE=R/4.5
SEEDS=(42 43 44 45 46)
WITHIN_ALPHA=1.0
BETWEEN_ALPHA=1.0
# =========================================================================

idx=$(( SLURM_ARRAY_TASK_ID - 1 ))
if [ "$idx" -lt 0 ] || [ "$idx" -ge "${#SEEDS[@]}" ]; then
    echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} has no entry in SEEDS (${#SEEDS[@]} seeds)." >&2
    exit 1
fi
export SMFSIM_SEED="${SEEDS[$idx]}"
BASE="${SMFSIM_SEEDS_OUTPUT:-results/bias_experiment_seeds}"
export SMFSIM_OUTPUT="${BASE}/seed_${SMFSIM_SEED}"

module load "$R_MODULE"
mkdir -p logs

echo "host          : $(hostname)"
echo "workdir       : $(pwd)"
echo "array task    : ${SLURM_ARRAY_TASK_ID} of ${#SEEDS[@]}"
echo "seed          : ${SMFSIM_SEED}"
echo "output        : ${SMFSIM_OUTPUT}"
echo "blocks        : ${SMFSIM_BLOCKS:-all}"
echo "extra args    : ${SMFSIM_EXTRA_ARGS:-none}"
echo "started       : $(date)"

echo "=== pre-flight ==="
Rscript -e '
  ok <- TRUE
  for (p in c("SMFnorm", "ComBatMet", "GenomicRanges", "data.table", "ggplot2")) {
    have <- requireNamespace(p, quietly = TRUE)
    cat(sprintf("  %-14s %s\n", p, if (have) "OK" else "MISSING"))
    if (!have) ok <- FALSE
  }
  if (!ok) quit(status = 1)
'

echo "=== running ==="
# shellcheck disable=SC2086  # SMFSIM_EXTRA_ARGS is word-split on purpose
Rscript SMFsim/inst/scripts/run_bias_experiment.R \
    --within_alpha  "$WITHIN_ALPHA" \
    --between_alpha "$BETWEEN_ALPHA" \
    ${SMFSIM_EXTRA_ARGS:-}

echo "finished      : $(date)"
