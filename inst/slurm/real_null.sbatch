#!/bin/bash
#SBATCH --job-name=smfsim_real_null
#SBATCH --account=jobrant          # <-- EDIT ME
#SBATCH --qos=jobrant-b            # <-- EDIT ME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jobrant@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=96gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Log paths are relative to the SUBMIT directory and are opened before this
# script runs: `mkdir -p logs` there first, and submit from the manuscript root.

# Real-data null: every 2-vs-2 split of one group's replicates, genome-wide,
# for raw, downsampled and several SMFnorm settings. Every DMR is a false
# positive. See SMFsim/inst/scripts/validate_real_null.R for outputs.
#
# Submit from the manuscript root:
#   sbatch SMFsim/inst/slurm/real_null.sbatch
# Options pass through the environment, e.g. chr1 only:
#   sbatch --export=ALL,REAL_NULL_CHR='^(chr)?1$' SMFsim/inst/slurm/real_null.sbatch
#
# Wall time is a guess: chr1 took ~10 minutes for 3 splits x 4 arms; the genome
# is ~11x the sites and this runs 6 arms.

set -euo pipefail

module load R/4.5

echo "host        : $(hostname)"
echo "workdir     : $(pwd)"
echo "group       : ${REAL_NULL_GROUP:-M1}"
echo "chromosomes : ${REAL_NULL_CHR:-all standard}"
echo "started     : $(date)"

Rscript SMFsim/inst/scripts/validate_real_null.R

echo "finished    : $(date)"
