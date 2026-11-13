#!/bin/bash
#SBATCH --job-name=smfsim_cbm_labels
#SBATCH --account=jobrant          # <-- EDIT ME
#SBATCH --qos=jobrant-b            # <-- EDIT ME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shenxinyu@ufl.edu,jobrant@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Log paths are relative to the SUBMIT directory and are opened before this
# script runs: submit from the manuscript root, where logs/ already exists.

# Is ComBatMet's null false-positive count a ComBat_met problem or an artifact
# of SMFsim's median-split batch labels? Reruns the matched-mean nulls under
# three labellings. See SMFsim/inst/scripts/diagnose_combatmet_batch_labels.R.
#
# Submit from the manuscript root:
#   sbatch SMFsim/inst/slurm/combatmet_batch_labels.sbatch
# Options go on the sbatch line, e.g. severe only:
#   sbatch --export=ALL,CBM_SCENARIOS=severe SMFsim/inst/slurm/combatmet_batch_labels.sbatch
#
# Time: ~15 min per ComBatMet run genome-wide; the default 3 scenarios x 3
# labellings is ~2.5 h.

set -euo pipefail

module load R/4.5

echo "host        : $(hostname)"
echo "workdir     : $(pwd)"
echo "scenarios   : ${CBM_SCENARIOS:-mild,moderate,severe}"
echo "chromosomes : ${CBM_CHR:-all standard}"
echo "output      : ${CBM_OUT:-results/combatmet_batch_labels}"
echo "started     : $(date)"

Rscript SMFsim/inst/scripts/diagnose_combatmet_batch_labels.R

echo "finished    : $(date)"
