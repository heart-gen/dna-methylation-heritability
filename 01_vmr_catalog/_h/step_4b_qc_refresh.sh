#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=vmr_qc_refresh
#SBATCH --output=logs/vmr_qc_refresh.%A.log
#
# 01_vmr_catalog step 4b: refresh QC tables on an already-accepted catalog.
#
# The accepted vmrcat-*-20260816 runs were sealed before any array probe
# universe existed on disk, so 04_turnover.R took its silent skip branch and
# qc/array_coverage.tsv held a placeholder note instead of the off-array numbers
# AGENTS.md 2.2 and 11 (Figure 1) require.
#
# Runs are immutable (AGENTS.md 5.2), so this mints a NEW run that copies the
# accepted catalog verbatim and re-runs QC only. The VMR calls are untouched and
# vmr_set_id carries forward unchanged, so nothing downstream is invalidated.
#
# Requires the 450K universe; build it first with
#   inputs/supportfiles/_h/step_1_build_array_universe.sh
#
# Usage, from the module's _m directory:
#   cd 01_vmr_catalog/_m && mkdir -p logs
#   COHORT=AA REGION=caudate SOURCE_RUN_ID=vmrcat-AA-caudate-20260816 \
#     sbatch ../_h/step_4b_qc_refresh.sh

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT=AA|all_individuals}"
: "${REGION:?set REGION}"
: "${SOURCE_RUN_ID:?set SOURCE_RUN_ID=<accepted vmrcat run>}"

require_file "$REPO_DIR/01_vmr_catalog/_m/runs/$SOURCE_RUN_ID/vmr/vmr_catalog.tsv"
require_file "$REPO_DIR/inputs/supportfiles/_m/450k_universe_hg38.tsv.gz"

log_job_info
log_message "**** Refreshing QC for ${COHORT}/${REGION} from ${SOURCE_RUN_ID} ****"

# 04b mints the run, copies the catalog tables, then invokes 04_turnover.R
# (array coverage) and 04c_genomic_context.R, and seals the result.
run_r "$REPO_DIR/01_vmr_catalog/_h/04b_rerun_array_coverage.R" \
    --cohort "$COHORT" --region "$REGION" --source-run-id "$SOURCE_RUN_ID"

log_message "**** Job ends ****"
