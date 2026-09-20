#!/bin/bash
# Submit one region of 09b_aging_application.
#
#   ./submit_aging.sh AA <region>
#
# Environment:
#   SMOKE_N=1   smoke run (unlocked config permitted, smoke permutation count)
#   AGE_RUN_ID_OVERRIDE=<id>   smoke only: choose the run ID
#   DRY_RUN=1   build the job graph, submit nothing
#
# Chain: open run (submit host) -> per-VMR age effects -> permuted axis test
#        -> coverage gate + region reading -> seal.
#
# No array: the age model is one matrix fit per spec over every VMR, so the
# region is one job. The cross-region stage (05) is run by hand after the three
# regions are accepted:
#   Rscript _h/05_cross_region_concordance.R --cohort AA

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/09b_aging_application"

COHORT=${1:?usage: submit_aging.sh AA <region>}
REGION=${2:?usage: submit_aging.sh AA <region>}

OPEN_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    OPEN_ARGS="--allow-unlocked"
    [ -n "${AGE_RUN_ID_OVERRIDE:-}" ] && OPEN_ARGS="$OPEN_ARGS --run-id ${AGE_RUN_ID_OVERRIDE}"
    log_message "SMOKE RUN: unlocked config permitted"
fi

RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
    --cohort "$COHORT" --region "$REGION" ${OPEN_ARGS} | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 00_new_run.R produced no run ID (an upstream gate most likely refused)" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
log_message "run ${RUN_ID}"

# Snapshot code and config, so an edit or branch switch cannot change a run in
# flight (see run_config.R).
mkdir -p "${RUN_DIR}/code"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"
mkdir -p "${RUN_DIR}/code/config"
cp -a "${REPO_DIR}/config/." "${RUN_DIR}/code/config/"
RUN_CODE="${RUN_DIR}/code/_h"

if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1, not submitting. Run dir: ${RUN_DIR}"
    cat <<GRAPH
planned job graph for ${RUN_ID}:
  1  01_age_effects.R    per-VMR age effects, every spec; checkpoint
  2  02_axis_test.R      afterok:1  permuted axis tests
  3  03_apply_gates.R    afterok:2
  4  04_finalize_run.R   afterok:3
GRAPH
    exit 0
fi

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
BASE_EXPORT="ALL,V2_REPO_ROOT=${REPO_DIR},V2_RUN_CODE=${RUN_CODE}"

sbatch_step () {  # name deps cpus mem time command
    local dep=()
    [ -n "$2" ] && dep=(--dependency="$2")
    sbatch --parsable "${dep[@]}" --chdir="${RUN_DIR}/logs" \
        --account=b1042 --partition=genomics --qos=buyin --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%j.out --error=%x-%j.err --export="$BASE_EXPORT" \
        --wrap="$6"
}
ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"

FIT_JOB=$(sbatch_step age-fit "" 2 16G 02:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/01_age_effects.R --run-id ${RUN_ID}")
printf '1\t01_age_effects.R\t%s\n' "$FIT_JOB" >> "$JOBS_TSV"

AXIS_JOB=$(sbatch_step age-axis "afterok:${FIT_JOB}" 2 16G 08:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/02_axis_test.R --run-id ${RUN_ID}")
printf '2\t02_axis_test.R\t%s\n' "$AXIS_JOB" >> "$JOBS_TSV"

GATE_JOB=$(sbatch_step age-gate "afterok:${AXIS_JOB}" 1 8G 00:30:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/03_apply_gates.R --run-id ${RUN_ID}")
printf '3\t03_apply_gates.R\t%s\n' "$GATE_JOB" >> "$JOBS_TSV"

FINAL_JOB=$(sbatch_step age-final "afterok:${GATE_JOB}" 1 8G 00:30:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/04_finalize_run.R --run-id ${RUN_ID}")
printf '4\t04_finalize_run.R\t%s\n' "$FINAL_JOB" >> "$JOBS_TSV"

log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
