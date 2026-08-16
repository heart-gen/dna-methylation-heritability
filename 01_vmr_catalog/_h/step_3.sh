#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=normal
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=vmr_summarize
#SBATCH --output=logs/vmr_summarize.%A.log
#
# 01_vmr_catalog step 3: call VMRs across all chromosomes, compute per-VMR
# methylation phenotypes, and mint the vmr_set_id. Not an array -- it must see
# every chromosome at once to reconcile them.

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT}"
: "${REGION:?set REGION}"
: "${RUN_ID:?set RUN_ID}"

log_job_info
log_message "**** Job starts: ${COHORT}/${REGION} run ${RUN_ID} ****"

run_r "$REPO_DIR/01_vmr_catalog/_h/02_summarize.R" \
    --cohort "$COHORT" --region "$REGION" --run-id "$RUN_ID" \
    ${ALLOW_UNLOCKED:+--allow-unlocked}

run_r "$REPO_DIR/01_vmr_catalog/_h/03_plot.R" \
    --cohort "$COHORT" --region "$REGION" --run-id "$RUN_ID" \
    ${ALLOW_UNLOCKED:+--allow-unlocked}

run_r "$REPO_DIR/01_vmr_catalog/_h/04_turnover.R" \
    --cohort "$COHORT" --region "$REGION" --run-id "$RUN_ID" \
    ${ALLOW_UNLOCKED:+--allow-unlocked}

log_message "**** Job ends ****"
