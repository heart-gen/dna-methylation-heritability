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

# SLURM copies the batch script to /var/spool, so ${BASH_SOURCE[0]} does NOT
# resolve to _h/ at run time. V2_REPO_ROOT is exported by the submit driver;
# fall back to the submit directory for a hand-run sbatch.
V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COHORT=${1:?usage: run_repeat_architecture.sh <AA|all_individuals>}
REGIONS=(caudate dlpfc hippocampus)

EXTRA_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    EXTRA_ARGS="--allow-unlocked"
    # SMOKE_N reaches the feature builder, which subsamples the catalog evenly.
    RRA_EXTRA_ARGS="--smoke-n ${SMOKE_N}"
    export RRA_EXTRA_ARGS
    log_message "SMOKE RUN: ${SMOKE_N} VMRs per cell, unlocked keys permitted"
fi

RUN_IDS=()
for REGION in "${REGIONS[@]}"; do
    log_message "opening ${COHORT} x ${REGION}"
    OPEN_ARGS="$EXTRA_ARGS"
    if [ -n "${SMOKE_N:-}" ] && [ -n "${RRA_RUN_ID_PREFIX:-}" ]; then
        OPEN_ARGS="$OPEN_ARGS --run-id ${RRA_RUN_ID_PREFIX}-${COHORT}-${REGION}-$(date +%Y%m%d)"
    fi
    RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
        --cohort "$COHORT" --region "$REGION" $OPEN_ARGS | tail -n 1)
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

# Feature building is a 64G, multi-hour job per cell -- it loads a BSseq object
# per chromosome for WGBS coverage. It runs on compute nodes, not here, and the
# three cells run in PARALLEL. Only the gate stage is cross-region, so only it
# waits for all three.
IFS=,; RUN_LIST="${RUN_IDS[*]}"; unset IFS
JOBS_TSV="${REPO_DIR}/04_repeat_repressive_architecture/_m/runs/${RUN_IDS[0]}/submitted-jobs.tsv"
printf 'step\tregion\tscript\tjob_id\n' > "$JOBS_TSV"

ASSOC_JOBS=()
for i in "${!RUN_IDS[@]}"; do
    RUN_ID="${RUN_IDS[$i]}"
    REGION="${REGIONS[$i]}"
    RUN_DIR="${REPO_DIR}/04_repeat_repressive_architecture/_m/runs/${RUN_ID}"
    EXPORT="ALL,RRA_RUN_ID=${RUN_ID},RRA_EXTRA_ARGS=${RRA_EXTRA_ARGS:-}"

    J1=$(sbatch --parsable --chdir="${RUN_DIR}/logs" --export="$EXPORT" \
        "${SCRIPT_DIR}/step_1_features.sh")
    printf '1\t%s\tstep_1_features.sh\t%s\n' "$REGION" "$J1" >> "$JOBS_TSV"

    J2=$(sbatch --parsable --dependency="afterok:${J1}" \
        --chdir="${RUN_DIR}/logs" --export="$EXPORT" \
        "${SCRIPT_DIR}/step_2_associate.sh")
    printf '2\t%s\tstep_2_associate.sh\t%s\n' "$REGION" "$J2" >> "$JOBS_TSV"
    ASSOC_JOBS+=("$J2")
done

IFS=:; DEP="${ASSOC_JOBS[*]}"; unset IFS
J3=$(sbatch --parsable --dependency="afterok:${DEP}" \
    --chdir="${REPO_DIR}/04_repeat_repressive_architecture/_m/runs/${RUN_IDS[0]}/logs" \
    --export="ALL,RRA_COHORT=${COHORT},RRA_RUN_IDS=${RUN_LIST}" \
    "${SCRIPT_DIR}/step_3_gates.sh")
printf '3\tall\tstep_3_gates.sh\t%s\n' "$J3" >> "$JOBS_TSV"

# Figures and sealing also span all three cells, and run only if the gates
# succeeded: a sealed cell whose claims table was never written would look
# finished while carrying no statement of what it may conclude.
J4=$(sbatch --parsable --dependency="afterok:${J3}" \
    --chdir="${REPO_DIR}/04_repeat_repressive_architecture/_m/runs/${RUN_IDS[0]}/logs" \
    --export="ALL,RRA_COHORT=${COHORT},RRA_RUN_IDS=${RUN_LIST}" \
    "${SCRIPT_DIR}/step_4_plot_finalize.sh")
printf '4\tall\tstep_4_plot_finalize.sh\t%s\n' "$J4" >> "$JOBS_TSV"

log_message "submitted 3 cells in parallel; gates ${J3} waits for all; plot/seal ${J4}"
log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
