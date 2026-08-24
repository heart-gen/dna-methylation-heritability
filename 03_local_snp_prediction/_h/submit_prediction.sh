#!/bin/bash
# Submit one cohort x region cell of 03_local_snp_prediction.
#
#   ./submit_prediction.sh <AA|all_individuals> <caudate|dlpfc|hippocampus>
#
# Environment:
#   VMRS_PER_ARRAY_TASK   VMRs per array task (default 5; nested CV is heavy)
#   MAX_CONCURRENT        array throttle (default 50, as 02 settled on)
#   SMOKE_N               run only the first N VMRs, with --allow-unlocked
#   DRY_RUN=1             build everything, submit nothing
#
# This driver refuses to submit unless 02 has an accepted run for the cell.
# Module 02 was rebuilt as the relative local-genetic-control pipeline on
# 2026-08-21 and its accepted-runs table is still empty, so the gate is CLOSED
# and this script has never been run in production.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/03_local_snp_prediction"

COHORT=${1:?usage: submit_prediction.sh <AA|all_individuals> <region>}
REGION=${2:?usage: submit_prediction.sh <AA|all_individuals> <region>}
VMRS_PER_ARRAY_TASK=${VMRS_PER_ARRAY_TASK:-5}
MAX_CONCURRENT=${MAX_CONCURRENT:-50}

EXTRA_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    # A smoke run is allowed to proceed on unlocked PI keys and an unaccepted
    # upstream, and is stamped smoke_run=TRUE in its manifest so it can never be
    # mistaken for production (AGENTS.md 14).
    EXTRA_ARGS="--allow-unlocked"
    log_message "SMOKE RUN: first ${SMOKE_N} VMRs, unlocked keys permitted"
fi

# ------------------------------------------------------------------ open run
log_message "opening run for ${COHORT} x ${REGION}"
RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
    --cohort "$COHORT" --region "$REGION" $EXTRA_ARGS | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 00_new_run.R produced no run ID (the upstream gate most likely refused)" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
log_message "run ${RUN_ID}"

# Snapshot the code into the run, so later edits to _h/ cannot change what an
# in-flight or already-finished run actually executed.
mkdir -p "${RUN_DIR}/code"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"

# --------------------------------------------------------------------- folds
run_r "${SCRIPT_DIR}/01_prepare_folds.R" --run-id "$RUN_ID" $EXTRA_ARGS

# ------------------------------------------------------------ chunk manifest
N_TASKS=$(( $(wc -l < "${RUN_DIR}/task-manifest.tsv") - 1 ))
if [ -n "${SMOKE_N:-}" ] && [ "$SMOKE_N" -lt "$N_TASKS" ]; then
    N_TASKS=$SMOKE_N
fi
N_CHUNKS=$(( (N_TASKS + VMRS_PER_ARRAY_TASK - 1) / VMRS_PER_ARRAY_TASK ))
CHUNK_MANIFEST="${RUN_DIR}/chunk-manifest.tsv"

awk -v n="$N_TASKS" -v size="$VMRS_PER_ARRAY_TASK" 'BEGIN {
    OFS = "\t"; print "chunk_id", "task_id"
    for (task = 1; task <= n; task++) print int((task - 1) / size) + 1, task
}' > "$CHUNK_MANIFEST"

# Validate what we just generated. A chunk manifest that silently drops or
# duplicates a task produces a run that reconciles as complete while missing
# VMRs, which is the failure mode reconcile() cannot catch on its own.
awk -F'\t' -v n="$N_TASKS" -v chunks="$N_CHUNKS" '
    NR == 1 { if ($1 != "chunk_id" || $2 != "task_id") exit 2; next }
    { rows++; seen[$2]++; if ($1 < 1 || $1 > chunks) exit 3 }
    END { if (rows != n) exit 4
          for (i = 1; i <= n; i++) if (seen[i] != 1) exit 5 }
' "$CHUNK_MANIFEST" || {
    echo "ERROR: generated chunk manifest failed validation" >&2; exit 1; }

printf 'expected_tasks\texpected_chunks\tvmrs_per_array_task\n%s\t%s\t%s\n' \
    "$N_TASKS" "$N_CHUNKS" "$VMRS_PER_ARRAY_TASK" > "${RUN_DIR}/expected-tasks.tsv"

log_message "${N_TASKS} VMRs in ${N_CHUNKS} chunks of ${VMRS_PER_ARRAY_TASK}"

if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1, not submitting. Run dir: ${RUN_DIR}"
    exit 0
fi

# -------------------------------------------------------------------- submit
JOB_ID=$(sbatch --parsable \
    --array="1-${N_CHUNKS}%${MAX_CONCURRENT}" \
    --chdir="${RUN_DIR}/logs" \
    --export=ALL,LSP_RUN_ID="$RUN_ID",LSP_CHUNK_MANIFEST="$CHUNK_MANIFEST",LSP_EXTRA_ARGS="$EXTRA_ARGS" \
    "${SCRIPT_DIR}/step_2_fit_oof.sh")

log_message "submitted array ${JOB_ID} (1-${N_CHUNKS}%${MAX_CONCURRENT})"
log_message "when it finishes: Rscript _h/03_combine_oof.R --run-id ${RUN_ID}"
