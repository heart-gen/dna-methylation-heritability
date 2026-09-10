#!/bin/bash
# Submit one cohort x region cell of 10_environmental_exploratory.
#
#   ./submit_environmental.sh <AA|all_individuals> <region>
#
# Environment:
#   SMOKE_N=1   smoke run (unlocked keys, unaccepted upstreams permitted)
#   DRY_RUN=1   build the job graph, submit nothing
#
# Chain: open run -> build and gate the exposure matrix -> [22 association
#        tasks] -> pool and correct once -> control-axis test -> coverage gate
#        -> seal.
#
# The exposure matrix is built on the submit host before anything is queued, on
# purpose: if no exposure clears the eligibility gate in this cell there is
# nothing to test, and that should cost seconds rather than a 22-task array.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/10_environmental_exploratory"

COHORT=${1:?usage: submit_environmental.sh <AA|all_individuals> <region>}
REGION=${2:?usage: submit_environmental.sh <AA|all_individuals> <region>}

OPEN_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    OPEN_ARGS="--allow-unlocked"
    [ -n "${ENV_RUN_ID_OVERRIDE:-}" ] && OPEN_ARGS="$OPEN_ARGS --run-id ${ENV_RUN_ID_OVERRIDE}"
    log_message "SMOKE RUN: unaccepted upstreams and unlocked keys permitted"
fi

RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
    --cohort "$COHORT" --region "$REGION" ${OPEN_ARGS} | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 00_new_run.R produced no run ID (an upstream gate most likely refused)" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
log_message "run ${RUN_ID}"

mkdir -p "${RUN_DIR}/code"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"
# Snapshot config alongside the code. Stages used to re-read config/ from the
# live working tree, so a config edit -- or a branch switch that removed the
# file -- changed or broke a run already in flight while its manifest attested
# to the original checksum.
mkdir -p "${RUN_DIR}/code/config"
cp -a "${REPO_DIR}/config/." "${RUN_DIR}/code/config/"
RUN_CODE="${RUN_DIR}/code/_h"

if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1, not submitting. Run dir: ${RUN_DIR}"
    cat <<GRAPH
planned job graph for ${RUN_ID}:
  1  01_build_exposure_matrix.R   submit host (donor set, composites, eligibility gate)
  2  step_1_association.sh        array 1-22: per-VMR exposure association
  3  02b_combine_associations.R   afterok:2  (one BH family per exposure, applied once)
  4  03_control_axis_test.R       afterok:3
  5  04_apply_gates.R             afterok:4
  6  05_finalize_run.R            afterok:5
GRAPH
    exit 0
fi

# The exposure matrix and its gate run on the submit host: seconds of work, and
# every array task reads the result. A cell with no eligible exposure stops
# here.
log_message "building the exposure matrix and applying the eligibility gate"
run_r "${RUN_CODE}/01_build_exposure_matrix.R" --run-id "$RUN_ID"

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
BASE_EXPORT="ALL,ENV_RUN_ID=${RUN_ID},V2_REPO_ROOT=${REPO_DIR},V2_RUN_CODE=${RUN_CODE}"

sbatch_step () {  # name deps cpus mem time command
    sbatch --parsable --dependency="$2" --chdir="${RUN_DIR}/logs" \
        --account=b1042 --partition=genomics --qos=buyin --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%j.out --error=%x-%j.err --export="$BASE_EXPORT" \
        --wrap="$6"
}
ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"

ASSOC_JOB=$(sbatch --parsable --chdir="${RUN_DIR}/logs" --export="$BASE_EXPORT" \
    "${RUN_CODE}/step_1_association.sh")
printf '2\tstep_1_association.sh\t%s\n' "$ASSOC_JOB" >> "$JOBS_TSV"

# Each BH family spans all autosomes, so it is corrected once, here, after the
# array completes -- never inside a per-chromosome task.
POOL_JOB=$(sbatch_step env-pool "afterok:${ASSOC_JOB}" 2 48G 02:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/02b_combine_associations.R --run-id ${RUN_ID}")
printf '3\t02b_combine_associations.R\t%s\n' "$POOL_JOB" >> "$JOBS_TSV"

AXIS_JOB=$(sbatch_step env-axis "afterok:${POOL_JOB}" 2 32G 01:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/03_control_axis_test.R --run-id ${RUN_ID}")
printf '4\t03_control_axis_test.R\t%s\n' "$AXIS_JOB" >> "$JOBS_TSV"

GATE_JOB=$(sbatch_step env-gate "afterok:${AXIS_JOB}" 1 16G 00:30:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/04_apply_gates.R --run-id ${RUN_ID}")
printf '5\t04_apply_gates.R\t%s\n' "$GATE_JOB" >> "$JOBS_TSV"

FINAL_JOB=$(sbatch_step env-final "afterok:${GATE_JOB}" 1 8G 00:30:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/05_finalize_run.R --run-id ${RUN_ID}")
printf '6\t05_finalize_run.R\t%s\n' "$FINAL_JOB" >> "$JOBS_TSV"

log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
