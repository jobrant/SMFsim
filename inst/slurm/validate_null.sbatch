#!/bin/bash
#SBATCH --job-name=smfsim_null_cal
#SBATCH --account=jobrant         # <-- EDIT ME
#SBATCH --qos=jobrant-b              # <-- EDIT ME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shenxinyu@ufl.edu,jobrant@ufl.edu 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err

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

mkdir -p LOGS
set -euo pipefail

module load R/4.5               # <-- EDIT ME, e.g. R/4.4

# Optional, passed with  sbatch --export=ALL,SMFSIM_OUTPUT=...,SMFSIM_EXTRA_ARGS='...'
#   SMFSIM_OUTPUT      output folder (default results/null_calibration)
#   SMFSIM_EXTRA_ARGS  extra parse_args options, e.g. '--bio_region_sd 0.5'
echo "host        : $(hostname)"
echo "workdir     : $(pwd)"
echo "output      : ${SMFSIM_OUTPUT:-results/null_calibration}"
echo "extra args  : ${SMFSIM_EXTRA_ARGS:-none}"
echo "started     : $(date)"

# shellcheck disable=SC2086  # SMFSIM_EXTRA_ARGS is word-split on purpose
Rscript /orange/jobrant/yuki/smfnorm-sim-manuscript/SMFsim/inst/scripts/validate_null_calibration.R ${SMFSIM_EXTRA_ARGS:-}

echo "finished    : $(date)"
echo
echo "Read the VERDICT block at the end of the .err file."
echo "Checks 1-3 PASS -> the significance path is calibrated and every method ran."
