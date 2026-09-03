#!/bin/bash
# Submit one cohort x region cell of 07_transcription_splicing_coupling.
#
#   ./submit_transcription_splicing.sh <AA|all_individuals> <region>
#
# Environment:
#   SMOKE_N=1   smoke run (unlocked keys, unaccepted upstreams permitted)
#   DRY_RUN=1   build the job graph, submit nothing
#
# Chain: open run -> build links (submit host) -> [one association job per
#        modality] -> coupling tests -> gate -> plot -> seal.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/07_transcription_splicing_coupling"

COHORT=${1:?usage: submit_transcription_splicing.sh <AA|all_individuals> <region>}
REGION=${2:?usage: submit_transcription_splicing.sh <AA|all_individuals> <region>}

PY_ENV="/projects/p32505/opt/envs/genomics"
run_py () { conda run --no-capture-output -p "$PY_ENV" python "$@"; }

OPEN_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    OPEN_ARGS="--allow-unlocked"
    [ -n "${TSC_RUN_ID_OVERRIDE:-}" ] && OPEN_ARGS="$OPEN_ARGS --run-id ${TSC_RUN_ID_OVERRIDE}"
    log_message "SMOKE RUN: unaccepted upstreams and unlocked keys permitted"
fi

MODALITIES=$(run_py -c "
import yaml
c = yaml.safe_load(open('${REPO_DIR}/config/transcription_splicing.yml'))
print(' '.join(k for k, v in c['modalities'].items() if v.get('enabled')))")
log_message "modalities: ${MODALITIES}"

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
  1  01_build_feature_links.R    submit host (tested universe on accepted VMRs)
  2  step_2_associations.sh      one job per modality: ${MODALITIES}
  3  03_test_coupling.R          afterok:2  (three tests per modality)
  4  04_apply_gates.R            afterok:3
  5  05_plot.py                  afterok:4
  6  06_finalize_run.R           afterok:5
GRAPH
    exit 0
fi

# Links are built on the submit host: seconds of work, and every modality job
# reads the output.
log_message "building the tested universe"
run_r "${RUN_CODE}/01_build_feature_links.R" --run-id "$RUN_ID"

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
BASE_EXPORT="ALL,TSC_RUN_ID=${RUN_ID},V2_REPO_ROOT=${REPO_DIR},V2_RUN_CODE=${RUN_CODE}"

ASSOC_DEPS=""
for MOD in $MODALITIES; do
    J=$(sbatch --parsable --chdir="${RUN_DIR}/logs" \
        --export="${BASE_EXPORT},TSC_MODALITY=${MOD}" \
        --job-name="tsc-assoc-${MOD}" \
        "${RUN_CODE}/step_2_associations.sh")
    printf '2\tstep_2_associations.sh:%s\t%s\n' "$MOD" "$J" >> "$JOBS_TSV"
    ASSOC_DEPS="${ASSOC_DEPS}:${J}"
done

sbatch_step () {  # name deps cpus mem time command
    sbatch --parsable --dependency="$2" --chdir="${RUN_DIR}/logs" \
        --account=b1042 --partition=genomics --qos=buyin --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%j.out --error=%x-%j.err --export="$BASE_EXPORT" \
        --wrap="$6"
}
R_ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"
PY_RUN="conda run --no-capture-output -p ${PY_ENV} python"

TEST_JOB=$(sbatch_step tsc-tests "afterok${ASSOC_DEPS}" 2 32G 02:00:00 \
    "${R_ENV_SRC} && run_r ${RUN_CODE}/03_test_coupling.R --run-id ${RUN_ID}")
printf '3\t03_test_coupling.R\t%s\n' "$TEST_JOB" >> "$JOBS_TSV"

GATE_JOB=$(sbatch_step tsc-gate "afterok:${TEST_JOB}" 1 16G 01:00:00 \
    "${R_ENV_SRC} && run_r ${RUN_CODE}/04_apply_gates.R --run-id ${RUN_ID}")
printf '4\t04_apply_gates.R\t%s\n' "$GATE_JOB" >> "$JOBS_TSV"

PLOT_JOB=$(sbatch_step tsc-plot "afterok:${GATE_JOB}" 1 16G 01:00:00 \
    "${R_ENV_SRC} && ${PY_RUN} ${RUN_CODE}/05_plot.py --run-id ${RUN_ID}")
printf '5\t05_plot.py\t%s\n' "$PLOT_JOB" >> "$JOBS_TSV"

FINAL_JOB=$(sbatch_step tsc-final "afterok:${PLOT_JOB}" 1 8G 01:00:00 \
    "${R_ENV_SRC} && run_r ${RUN_CODE}/06_finalize_run.R --run-id ${RUN_ID}")
printf '6\t06_finalize_run.R\t%s\n' "$FINAL_JOB" >> "$JOBS_TSV"

log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
