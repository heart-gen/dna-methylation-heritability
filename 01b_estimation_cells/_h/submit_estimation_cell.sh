#!/bin/bash
# Materialize one donor-group estimation cell, end to end.
#
#   cd 01b_estimation_cells/_m && mkdir -p logs
#   COHORT=all_individuals GROUP=EA REGION=dlpfc \
#     VMR_RUN_ID=vmrcat-all_individuals-dlpfc-20260816 \
#     ../_h/submit_estimation_cell.sh
#
# Stage 00 runs on the SUBMIT HOST, not under sbatch: it is what mints the run
# ID that every later stage needs, and the VMR count it copies is what sizes the
# stage-02 array. Stages 01-04 are chained with SLURM dependencies.
#
# Stage 04 depends with afterok, not afterany: unlike a reconciliation stage,
# closing a run that never finished extracting would seal a partial cell.
# Stage 02's own reconciliation is inside stage 04, and a cancelled array leaves
# VMRs with neither a BED nor a marker, which stage 04 reports and refuses.
#
# Set DRY_RUN=TRUE to print the plan without submitting.

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT=all_individuals}"
: "${GROUP:?set GROUP=AA|EA}"
: "${REGION:?set REGION}"
: "${VMR_RUN_ID:?set VMR_RUN_ID to an accepted 01_vmr_catalog run}"

H_DIR="$REPO_DIR/01b_estimation_cells/_h"
DRY_RUN="${DRY_RUN:-FALSE}"
THROTTLE="${THROTTLE:-250}"
OPEN_ARGS=""
if [ "${ALLOW_UNLOCKED:-FALSE}" = "TRUE" ]; then
    OPEN_ARGS="--allow-unlocked"
fi

echo "cell      ${COHORT}.${GROUP} / ${REGION}"
echo "catalog   ${VMR_RUN_ID}"

if [ "$DRY_RUN" = "TRUE" ]; then
    echo "DRY_RUN: would open the run, then chain steps 1-4."
    exit 0
fi

# Stage 00 prints the run ID as its last line.
RUN_ID=$(run_r "$H_DIR/00_new_run.R" \
    --cohort "$COHORT" --group "$GROUP" --region "$REGION" \
    --vmr-run-id "$VMR_RUN_ID" $OPEN_ARGS | tail -n 1)
: "${RUN_ID:?stage 00 produced no run ID}"

RUN_DIR="$REPO_DIR/01b_estimation_cells/_m/runs/$RUN_ID"
require_file "$RUN_DIR/manifest.tsv"
echo "run       $RUN_ID"

N=$(wc -l < "$RUN_DIR/vmr/vmr.bed")
echo "VMRs      $N (throttle %${THROTTLE})"

export RUN_ID V2_REPO_ROOT

JOB_PCA=$(sbatch --parsable "$H_DIR/step_1_group_pca.sh")
JOB_EXTRACT=$(sbatch --parsable --array=1-${N}%${THROTTLE} \
    "$H_DIR/step_2_extract_cis.sh")
JOB_COVS=$(sbatch --parsable --dependency=afterok:"$JOB_EXTRACT" \
    "$H_DIR/step_3_subset_covariates.sh")
JOB_CLOSE=$(sbatch --parsable \
    --dependency=afterok:"$JOB_PCA":"$JOB_COVS" "$H_DIR/step_4_close.sh")

{
    printf 'step\tjob_id\n'
    printf 'step_1_group_pca\t%s\n' "$JOB_PCA"
    printf 'step_2_extract_cis\t%s\n' "$JOB_EXTRACT"
    printf 'step_3_subset_covariates\t%s\n' "$JOB_COVS"
    printf 'step_4_close\t%s\n' "$JOB_CLOSE"
} > "$RUN_DIR/submitted-jobs.tsv"

cat "$RUN_DIR/submitted-jobs.tsv"
echo "Nothing is accepted until step 4 seals the run AND a human records it in"
echo "01b_estimation_cells/README.md (AGENTS.md 6)."
