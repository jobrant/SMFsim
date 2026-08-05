#!/bin/bash
#SBATCH --job-name=smfsim_bias
#SBATCH --account=EDIT_ME_ACCOUNT          # <-- EDIT ME
#SBATCH --qos=EDIT_ME_QOS                  # <-- EDIT ME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jobrant@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=96:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# NOTE: log paths are relative to the SUBMIT directory, not $WORKDIR, and SLURM
# opens them before this script runs. Run `mkdir -p logs` in your submit
# directory first, and submit from the same place every time.

# THE HEADLINE RUN. Produces Figures 1 and 2.
#
# Blocks (all genome-wide, all four methods):
#   control_within   specificity control - matched means, no artifact
#   bias_parametric  the artifact - systematic between-group efficiency bias
#
# PREREQUISITES, in order:
#   1. The alpha grid has been stitched and a pair LOCKED. Set them below.
#   2. validate_null_calibration.R passed genome-wide (FAST_MODE <- FALSE).
#   3. ComBatMet is installed (checked below - this run is the first time it
#      contributes real numbers, so a missing package wastes the whole job).
#   4. Any previous results/bias_experiment has been archived. The script's
#      overwrite guard will abort rather than clobber it:
#        mv results/bias_experiment results/bias_experiment_$(date +%F)_superseded
#
# Wall time is a guess: the genome-wide alpha sweep took roughly 8 h per cell
# with ONE method; this runs four. Raise --time or split the blocks if it is
# tight, and check your QOS ceiling.

set -euo pipefail

# === EDIT ME =============================================================
WORKDIR=EDIT_ME_PATH_TO_SMFsim             # package root
R_MODULE=EDIT_ME_R_MODULE                  # e.g. R/4.4

# The LOCKED alphas from the stitched grid. Leave as UNSET to abort.
WITHIN_ALPHA=UNSET
BETWEEN_ALPHA=UNSET
# =========================================================================

# --- Refuse to burn 4 days at the wrong parameters -----------------------
if [ "$WITHIN_ALPHA" = "UNSET" ] || [ "$BETWEEN_ALPHA" = "UNSET" ]; then
    echo "ERROR: WITHIN_ALPHA / BETWEEN_ALPHA are still UNSET." >&2
    echo "Set them to the values locked from the stitched alpha grid" >&2
    echo "(see the 'Separability of the two alphas' block in the stitch log)." >&2
    echo >&2
    echo "Running at the parse_args defaults would silently produce results at" >&2
    echo "parameters you did not choose - which is exactly how the June figures" >&2
    echo "ended up superseded." >&2
    exit 1
fi

module load "$R_MODULE"
cd "$WORKDIR"
mkdir -p logs

echo "host          : $(hostname)"
echo "workdir       : $(pwd)"
echo "within_alpha  : ${WITHIN_ALPHA}"
echo "between_alpha : ${BETWEEN_ALPHA}"
echo "started       : $(date)"
echo

# --- Pre-flight: fail in seconds, not on day three -----------------------
echo "=== pre-flight ==="
Rscript -e '
  ok <- TRUE
  for (p in c("SMFnorm", "ComBatMet", "data.table", "ggplot2")) {
    have <- requireNamespace(p, quietly = TRUE)
    cat(sprintf("  %-12s %s\n", p, if (have) "OK" else "MISSING"))
    if (!have) ok <- FALSE
  }
  if (!ok) {
    cat("\nInstall the missing package(s) before submitting.\n")
    cat("  ComBatMet: remotes::install_github(\"JmWangBio/ComBatMet\")\n")
    quit(status = 1)
  }
'

if [ -d results/bias_experiment ] && [ -n "$(ls -A results/bias_experiment)" ]; then
    echo
    echo "ERROR: results/bias_experiment already exists and is not empty." >&2
    echo "The script's guard will abort. Archive it first:" >&2
    echo "  mv results/bias_experiment results/bias_experiment_\$(date +%F)_superseded" >&2
    exit 1
fi

echo
echo "=== running ==="
Rscript inst/scripts/run_bias_experiment.R \
    --within_alpha  "$WITHIN_ALPHA" \
    --between_alpha "$BETWEEN_ALPHA"

echo
echo "finished      : $(date)"
cat <<'EOF'

Outputs under results/bias_experiment/:
  control_within/   specificity control
  bias_parametric/  the artifact -> Figures 1 and 2
  each block: null_simulation_results.csv, spikein_simulation_results.csv,
              <block>_results.rds (holds the EFFECTIVE config), figures/

ARCHIVE THIS RUN before anything else writes to that path:
  cp -r results/bias_experiment results/bias_experiment_$(date +%F)_locked

Confirm the run used the intended parameters before trusting the figures:
  r <- readRDS("results/bias_experiment/bias_parametric/bias_parametric_results.rds")
  r$config[c("within_alpha", "between_alpha", "min_coverage", "dispersion_s",
             "metilene_q_source", "metilene_p_column", "metilene_mtc")]
EOF
