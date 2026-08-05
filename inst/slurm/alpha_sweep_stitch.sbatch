#!/bin/bash
#SBATCH --job-name=smfsim_alpha_stitch
#SBATCH --account=EDIT_ME_ACCOUNT          # <-- EDIT ME
#SBATCH --qos=EDIT_ME_QOS                  # <-- EDIT ME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jobrant@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# NOTE: relative to the SUBMIT directory, not $WORKDIR, and opened before this
# script runs -- `mkdir -p logs` below is too late to help. Create it in your
# submit directory first. The separability verdict goes to stdout, so losing
# this file means losing the main result of the stitch step.

# Stitch the array results together, then report and plot.
#
# Runs run_alpha_sweep.R with NO SLURM_ARRAY_TASK_ID set, so it walks the whole
# grid, finds every cell already complete, reads them from cache, and produces:
#   alpha_sweep_null.csv / alpha_sweep_spikein.csv   (combined)
#   alpha_grid_F1.png / alpha_grid_null_FP.png       (heatmaps)
#   a SEPARABILITY report on stdout
#
# Submit with a dependency so it only runs if every cell succeeded:
#   ARRAY_JOB=$(sbatch --parsable inst/slurm/alpha_sweep_array.sbatch)
#   sbatch --dependency=afterok:$ARRAY_JOB inst/slurm/alpha_sweep_stitch.sbatch
#
# If some cells failed, fix and resubmit the array first -- completed cells are
# skipped on resume, so only the missing ones re-run. Running this job while
# cells are missing silently produces a partial grid.

set -euo pipefail

WORKDIR=EDIT_ME_PATH_TO_SMFsim             # <-- EDIT ME (package root)
module load EDIT_ME_R_MODULE               # <-- EDIT ME, e.g. R/4.4

cd "$WORKDIR"
mkdir -p logs

# Guard against stitching an incomplete grid. Expected cell count must match
# the --array range used for the sweep.
EXPECTED_CELLS=15
SWEEP_DIR=results/alpha_sweep_fast          # results/alpha_sweep if FAST_TUNE=FALSE

if [ -d "$SWEEP_DIR" ]; then
    FOUND=$(find "$SWEEP_DIR" -mindepth 2 -maxdepth 2 \
                 -name 'spikein_simulation_results.csv' | wc -l)
    echo "completed cells: ${FOUND} of ${EXPECTED_CELLS}"
    if [ "$FOUND" -lt "$EXPECTED_CELLS" ]; then
        echo "WARNING: grid is incomplete. The stitched CSVs and the"
        echo "separability verdict will be based on a PARTIAL grid."
        echo "Resubmit the array to fill the gaps before trusting the result."
    fi
else
    echo "WARNING: ${SWEEP_DIR} not found - did the array job run?"
fi

echo "started     : $(date)"

Rscript inst/scripts/run_alpha_sweep.R

echo "finished    : $(date)"
echo
echo "Look for the 'Separability of the two alphas' block in this log:"
echo "  SEPARABLE     -> the two alphas can be tuned and reported independently"
echo "  NOT SEPARABLE -> quote the joint optimum, not two 1-D sweeps"
