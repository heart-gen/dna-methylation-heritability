#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=vmr_close_run
#SBATCH --output=logs/vmr_close_run.%A.log
#
# 01_vmr_catalog step 5: seal the run.
#
# MUST be the last step. It checksums every output and makes the run directory
# read-only, so it has to run after step_4 has written plink_format/. Closing
# earlier would leave step_4 unable to write into its own run.

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT}"
: "${REGION:?set REGION}"
: "${RUN_ID:?set RUN_ID}"

log_job_info
log_message "**** Sealing run ${RUN_ID} ****"

run_r "$REPO_DIR/01_vmr_catalog/_h/05_close_run.R" \
    --cohort "$COHORT" --region "$REGION" --run-id "$RUN_ID" \
    ${ALLOW_UNLOCKED:+--allow-unlocked}

log_message "**** Job ends ****"
