#!/bin/bash
# Run 04_repeat_repressive_architecture for one cohort, all three regions.
#
#   ./run_repeat_architecture.sh <AA|all_individuals>
#
# 04 is not an array workload -- it is three per-region model fits over a few
# thousand rows each -- so it runs as a single batch job rather than a SLURM
# array. All three regions run in one job because the interpretation gates
# (03_apply_gates.R) are defined ACROSS regions and cannot be applied to a
# partial set.
#
# Environment:
#   SMOKE_N=1   permit unlocked keys and unaccepted upstreams
#   DRY_RUN=1   build run directories, fit nothing

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COHORT=${1:?usage: run_repeat_architecture.sh <AA|all_individuals>}
REGIONS=(caudate dlpfc hippocampus)

EXTRA_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    EXTRA_ARGS="--allow-unlocked"
    log_message "SMOKE RUN: unaccepted upstreams and unlocked keys permitted"
fi

RUN_IDS=()
for REGION in "${REGIONS[@]}"; do
    log_message "opening ${COHORT} x ${REGION}"
    RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
        --cohort "$COHORT" --region "$REGION" $EXTRA_ARGS | tail -n 1)
    if [ -z "$RUN_ID" ]; then
        echo "ERROR: could not open a run for ${REGION} (upstream gate refused)" >&2
        exit 1
    fi
    RUN_DIR="${REPO_DIR}/04_repeat_repressive_architecture/_m/runs/${RUN_ID}"
    mkdir -p "${RUN_DIR}/code"
    cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"
    RUN_IDS+=("$RUN_ID")
done

if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1, stopping after run creation: ${RUN_IDS[*]}"
    exit 0
fi

for RUN_ID in "${RUN_IDS[@]}"; do
    log_message "features + models for ${RUN_ID}"
    run_r "${SCRIPT_DIR}/01_build_features.R"   --run-id "$RUN_ID"
    run_r "${SCRIPT_DIR}/02_test_association.R" --run-id "$RUN_ID"
done

IFS=,; RUN_LIST="${RUN_IDS[*]}"; unset IFS
log_message "applying cross-region interpretation gates"
run_r "${SCRIPT_DIR}/03_apply_gates.R" --cohort "$COHORT" --run-ids "$RUN_LIST"

log_message "done: ${RUN_LIST}"
