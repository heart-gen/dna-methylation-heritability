#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --array=0-1
#SBATCH --job-name=vmr_prepare_sex
#SBATCH --output=logs/vmr_prepare_sex.%A_%a.log
#
# 01_vmr_catalog step 1x: sex chromosomes, reported separately.
#
# AGENTS.md 7.1: "sex chromosomes reported separately or excluded with an
# explicit manifest". These CpGs have no C->T SNP mask, which is why the legacy
# caudate catalog carried a 3x excess of X/Y VMRs (431 vs 143/147) -- defect V4.
# Output goes to _m/runs/{RUN_ID}/excluded/ and never joins the primary catalog.
#
# Optional: run it only if you need the sex-chromosome manifest.

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT}"
: "${REGION:?set REGION}"
: "${RUN_ID:?set RUN_ID}"

SEX_CHROMS=(X Y)
CHR="${SEX_CHROMS[${SLURM_ARRAY_TASK_ID:-0}]}"

log_job_info
log_message "**** Sex chromosome chr${CHR} (excluded from primary catalog) ****"

run_r "$REPO_DIR/01_vmr_catalog/_h/00_prepare.R" \
    --cohort "$COHORT" --region "$REGION" --chrom "$CHR" --run-id "$RUN_ID" \
    ${ALLOW_UNLOCKED:+--allow-unlocked}

log_message "**** Job ends ****"
