#!/bin/bash
# Run 08_region_donor_generalization end to end.
#
#   cd 08_region_donor_generalization/_m && mkdir -p logs
#   COHORT=AA ../_h/submit_region_donor_generalization.sh
#
# Every stage is cheap -- this module assembles accepted upstream results and
# never refits a model -- so the whole chain runs on the submit host rather than
# through a SLURM array. The expensive part of tier 3 is the 02/03 runs on the
# AA.n118r* cells, which are separate modules' work and must be sealed and
# accepted before this runs at all.
#
# Set DRY_RUN=TRUE to print the plan without opening a run.

set -euo pipefail

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT=AA}"

H_DIR="$REPO_DIR/08_region_donor_generalization/_h"
DRY_RUN="${DRY_RUN:-FALSE}"
OPEN_ARGS=""
if [ "${ALLOW_UNLOCKED:-FALSE}" = "TRUE" ]; then
    OPEN_ARGS="--allow-unlocked"
fi

echo "cohort    $COHORT"
echo "stages    00 open, 01 tier1, 02 tier2+tier4, 03 donor group, 04 tier3, 05 gate, 06 seal"

if [ "$DRY_RUN" = "TRUE" ]; then
    echo "DRY_RUN: would open a run and chain stages 01-06."
    exit 0
fi

# Stage 00 prints the run ID as its last line.
RUN_ID=$(run_r "$H_DIR/00_new_run.R" --cohort "$COHORT" $OPEN_ARGS | tail -n 1)
: "${RUN_ID:?stage 00 produced no run ID}"

RUN_DIR="$REPO_DIR/08_region_donor_generalization/_m/runs/$RUN_ID"
require_file "$RUN_DIR/manifest.tsv"
echo "run       $RUN_ID"

for stage in 01_cross_region_replication 02_identified_difference \
             03_donor_group_concordance 04_caudate_downsampling \
             05_apply_gates 06_finalize_run; do
    echo "=== $stage ==="
    run_r "$H_DIR/${stage}.R" --run-id "$RUN_ID"
done

echo
echo "Sealed $RUN_ID. Nothing is accepted until a human records it in"
echo "08_region_donor_generalization/README.md (AGENTS.md 6)."
