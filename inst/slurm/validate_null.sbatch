#!/bin/bash
#SBATCH --job-name=smfsim_null_cal
#SBATCH --account=EDIT_ME_ACCOUNT          # <-- EDIT ME
#SBATCH --qos=EDIT_ME_QOS                  # <-- EDIT ME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jobrant@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48gb
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Calibration gate. Run this BEFORE the alpha grid.
#
# Checks, in the order they can invalidate downstream work:
#   1. metilene's "Number of Tests" line is being parsed (BH denominator)
#   2. a TRUE NULL yields ~0 false positives
#   3. the rewritten ComBatMet runs end-to-end
#
# Run FAST_MODE <- TRUE first (chr1, edit the script), confirm the verdict, then
# set FAST_MODE <- FALSE and resubmit for the genome-wide answer. The 12 h wall
# time is sized for the genome-wide run; chr1 finishes in well under an hour.

set -euo pipefail

WORKDIR=EDIT_ME_PATH_TO_SMFsim             # <-- EDIT ME (package root)
module load EDIT_ME_R_MODULE               # <-- EDIT ME, e.g. R/4.4

cd "$WORKDIR"
mkdir -p logs

echo "host        : $(hostname)"
echo "workdir     : $(pwd)"
echo "started     : $(date)"

Rscript inst/scripts/validate_null_calibration.R

echo "finished    : $(date)"
echo
echo "Read the VERDICT block at the end of the .out file."
echo "All three PASS -> proceed to alpha_sweep_array.sbatch."
echo "Any FAIL       -> stop and fix; do not launch the grid."
