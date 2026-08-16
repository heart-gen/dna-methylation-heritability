#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=normal
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# chr22 (489,722 CpGs x 153 donors) peaked at 8.3 GB. chr1 carries roughly 5x
# the CpGs and 01_analyze.R holds the full matrix in memory, so budget ~40 GB
# there. 96G is what the chr22 smoke run was allocated and verified at.
#SBATCH --mem=96G
#SBATCH --mail-type=FAIL
#SBATCH --array=1-22
#SBATCH --job-name=vmr_analyze
#SBATCH --output=logs/vmr_analyze.%A_%a.log
#
# 01_vmr_catalog step 2: methylation PCs and residual variance per chromosome.
# This is the step that carries the V1 alignment fix.

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT}"
: "${REGION:?set REGION}"
: "${RUN_ID:?set RUN_ID}"

CHR="${SLURM_ARRAY_TASK_ID:-${CHR:?set CHR when running outside SLURM}}"

log_job_info
log_message "**** Job starts: ${COHORT}/${REGION} chr${CHR} run ${RUN_ID} ****"

run_r "$REPO_DIR/01_vmr_catalog/_h/01_analyze.R" \
    --cohort "$COHORT" --region "$REGION" --chrom "$CHR" --run-id "$RUN_ID" \
    ${ALLOW_UNLOCKED:+--allow-unlocked}

log_message "**** Job ends ****"
