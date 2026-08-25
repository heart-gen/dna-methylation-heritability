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
#   LSP_RUN_ID_OVERRIDE   explicit run ID; smoke runs only
#
# This driver refuses to submit unless 02 has an accepted run for the cell.
# Module 02 recorded six accepted runs on 2026-08-23, so the gate is OPEN for
# every cohort x region cell.
#
# Stages chain by SLURM dependency: 1 -> 2 (array) -> 3 -> 4 -> 5. Step 3
# depends with afterany so a scheduler cancellation is RECONCILED rather than
# silently leaving the audit blocked; 4 and 5 use afterok, because there is no
# sense in gating or sealing a run whose combine step failed.

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
    OPEN_ARGS="--allow-unlocked --smoke-n ${SMOKE_N}"
    if [ -n "${LSP_RUN_ID_OVERRIDE:-}" ]; then
        OPEN_ARGS="$OPEN_ARGS --run-id ${LSP_RUN_ID_OVERRIDE}"
    fi
    log_message "SMOKE RUN: ${SMOKE_N} VMRs spread across the catalog, unlocked keys permitted"
fi

# ------------------------------------------------------------------ open run
log_message "opening run for ${COHORT} x ${REGION}"
OPEN_ARGS=${OPEN_ARGS:-$EXTRA_ARGS}
RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
    --cohort "$COHORT" --region "$REGION" $OPEN_ARGS | tail -n 1)
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

# Everything from here on executes the SNAPSHOT, not the live _h/. The copy
# above is only provenance if the jobs actually run it: previously the drivers
# snapshotted _h/ and then sbatch'd ${SCRIPT_DIR}, so an edit to _h/ while a run
# was queued changed what that run executed while its manifest and snapshot
# attested to the older commit.
RUN_CODE="${RUN_DIR}/code/_h"

# Folds are assigned by step 1 as a scheduled job, not inline here, so the
# whole pipeline is reproducible from the recorded job graph rather than partly
# from whatever the submit host happened to do.

# ------------------------------------------------------------ chunk manifest
# The task manifest is authoritative: 00_new_run.R --smoke-n has already
# shrunk it for a smoke run. Truncating again here would leave tasks expected
# but unchunked, and reconcile() would (correctly) fail the run.
N_TASKS=$(( $(wc -l < "${RUN_DIR}/task-manifest.tsv") - 1 ))
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
    cat <<GRAPH
planned job graph for ${RUN_ID}:
  1  step_1_prepare_folds.sh   -
  2  step_2_fit_oof.sh         afterok:1, array 1-${N_CHUNKS}%${MAX_CONCURRENT}
  3  step_3_combine_oof.sh     afterany:2   (reconciles cancelled array tasks)
  4  step_4_check.sh           afterok:3
  5  step_5_finalize.sh        afterok:4
GRAPH
    exit 0
fi

# -------------------------------------------------------------------- submit
JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
COMMON_EXPORT="ALL,LSP_RUN_ID=${RUN_ID},LSP_EXTRA_ARGS=${EXTRA_ARGS},V2_RUN_CODE=${RUN_CODE}"

JOB_FOLDS=$(sbatch --parsable --chdir="${RUN_DIR}/logs" \
    --export="$COMMON_EXPORT" "${RUN_CODE}/step_1_prepare_folds.sh")
printf '1\tstep_1_prepare_folds.sh\t%s\n' "$JOB_FOLDS" >> "$JOBS_TSV"

JOB_FIT=$(sbatch --parsable \
    --array="1-${N_CHUNKS}%${MAX_CONCURRENT}" \
    --dependency="afterok:${JOB_FOLDS}" \
    --chdir="${RUN_DIR}/logs" \
    --export="${COMMON_EXPORT},LSP_CHUNK_MANIFEST=${CHUNK_MANIFEST}" \
    "${RUN_CODE}/step_2_fit_oof.sh")
printf '2\tstep_2_fit_oof.sh\t%s\n' "$JOB_FIT" >> "$JOBS_TSV"

JOB_COMB=$(sbatch --parsable --dependency="afterany:${JOB_FIT}" \
    --chdir="${RUN_DIR}/logs" --export="$COMMON_EXPORT" \
    "${RUN_CODE}/step_3_combine_oof.sh")
printf '3\tstep_3_combine_oof.sh\t%s\n' "$JOB_COMB" >> "$JOBS_TSV"

JOB_CHECK=$(sbatch --parsable --dependency="afterok:${JOB_COMB}" \
    --chdir="${RUN_DIR}/logs" --export="$COMMON_EXPORT" \
    "${RUN_CODE}/step_4_check.sh")
printf '4\tstep_4_check.sh\t%s\n' "$JOB_CHECK" >> "$JOBS_TSV"

JOB_FINAL=$(sbatch --parsable --dependency="afterok:${JOB_CHECK}" \
    --chdir="${RUN_DIR}/logs" --export="$COMMON_EXPORT" \
    "${RUN_CODE}/step_5_finalize.sh")
printf '5\tstep_5_finalize.sh\t%s\n' "$JOB_FINAL" >> "$JOBS_TSV"

log_message "submitted array ${JOB_FIT} (1-${N_CHUNKS}%${MAX_CONCURRENT}) and the 1->5 chain"
log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
