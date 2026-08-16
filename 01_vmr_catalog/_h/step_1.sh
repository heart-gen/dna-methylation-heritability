#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --array=1-22
#SBATCH --job-name=vmr_prepare
#SBATCH --output=logs/vmr_prepare.%A_%a.log
#
# 01_vmr_catalog step 1: per-chromosome CpG matrix, covariates, donor manifest.
#
# Autosomes only (--array=1-22). Sex chromosomes are handled by step_1x.sh and
# land in the run's excluded/ directory -- defect V4.
#
# Submit from the module's _m/ directory:
#   cd 01_vmr_catalog/_m && mkdir -p logs
#   COHORT=AA REGION=caudate RUN_ID=<id> sbatch ../_h/step_1.sh

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT=AA|all_individuals}"
: "${REGION:?set REGION=caudate|dlpfc|hippocampus}"
: "${RUN_ID:?set RUN_ID (create it with 00_new_run.R)}"

CHR="${SLURM_ARRAY_TASK_ID:-${CHR:?set CHR when running outside SLURM}}"

log_job_info
log_message "**** Job starts: ${COHORT}/${REGION} chr${CHR} run ${RUN_ID} ****"

run_r "$REPO_DIR/01_vmr_catalog/_h/00_prepare.R" \
    --cohort "$COHORT" --region "$REGION" --chrom "$CHR" --run-id "$RUN_ID" \
    ${ALLOW_UNLOCKED:+--allow-unlocked}

log_message "**** Job ends ****"
