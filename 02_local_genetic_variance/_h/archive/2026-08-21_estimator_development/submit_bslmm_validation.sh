#!/bin/bash
#
# Full locked Module 02 validation grid for raw BSLMM PVE.
# Does NOT alter observed-data production. Replacement requires PASS on
# config/acceptance-criteria.tsv plus domain/failure gates.
#
# Usage, from module _m/:
#   ../_h/submit_bslmm_validation.sh [RUN_ID] [CONFIG_TSV]

set -euo pipefail

usage() {
    echo "Usage: $0 [RUN_ID] [CONFIG_TSV]"
    echo "Default config: config/bslmm-validation.tsv"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# > 2 )); then usage >&2; exit 1; fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
DRY_RUN=${SUBMIT_BSLMM_VAL_DRY_RUN:-FALSE}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
KEEP_WORK=${BSLMM_VAL_KEEP_WORK:-FALSE}

CONFIG=${2:-${ANALYSIS_DIR}/config/bslmm-validation.tsv}
CONFIG=$(readlink -f "${CONFIG}")
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

RUN_ID=${1:-bslmm-validation-$(date -u +%Y%m%d)}
[[ "${RUN_ID}" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Invalid RUN_ID" >&2; exit 1; }

[[ -x "${ENV_PATH}/bin/Rscript" ]] || { echo "Missing Rscript in ${ENV_PATH}" >&2; exit 1; }
export PATH="/projects/p32505/opt/bin:${PATH}"
command -v gemma >/dev/null || { echo "gemma not on PATH" >&2; exit 1; }
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

MAX_CONCURRENT=${MAX_CONCURRENT:-$(awk -F '\t' '$1=="max_concurrent"{print $2; exit}' "${CONFIG}")}
SCENARIOS_PER_ARRAY_TASK=${SCENARIOS_PER_ARRAY_TASK:-$(awk -F '\t' '$1=="scenarios_per_array_task"{print $2; exit}' "${CONFIG}")}

RUN_ROOT=${RUN_BASE}/${RUN_ID}
if [[ -e "${RUN_ROOT}" ]]; then
    echo "Run directory already exists: ${RUN_ROOT}" >&2
    exit 1
fi
mkdir -p "${RUN_ROOT}/config" "${RUN_ROOT}/logs" "${RUN_ROOT}/provenance" \
    "${RUN_ROOT}/code" "${RUN_ROOT}/raw" "${RUN_ROOT}/combined" "${RUN_ROOT}/work"

cp "${CONFIG}" "${RUN_ROOT}/config/bslmm-validation.tsv"
cp "${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
   "${RUN_ROOT}/config/acceptance-criteria.tsv"
cp -a "${SCRIPT_DIR}" "${RUN_ROOT}/code/_h"
# Refresh helpers from live tree if snapshot raced with edits
cp "${SCRIPT_DIR}/bslmm_pilot_functions.R" "${RUN_ROOT}/code/_h/bslmm_pilot_functions.R"
cp "${SCRIPT_DIR}/16_run_bslmm_validation_scenario.R" "${RUN_ROOT}/code/_h/"
cp "${SCRIPT_DIR}/17_evaluate_bslmm_validation.R" "${RUN_ROOT}/code/_h/"
cp "${SCRIPT_DIR}/15_make_bslmm_validation_manifest.R" "${RUN_ROOT}/code/_h/"
RUN_SCRIPT_DIR=${RUN_ROOT}/code/_h
MANIFEST=${RUN_ROOT}/config/scenarios.tsv

cd "${REPO_ROOT}"
"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/15_make_bslmm_validation_manifest.R" \
    --config="${RUN_ROOT}/config/bslmm-validation.tsv" \
    --output="${MANIFEST}"
TASKS=$(($(wc -l < "${MANIFEST}") - 1))
(( TASKS >= 1 )) || { echo "Empty manifest" >&2; exit 1; }

CHUNK_MANIFEST=${RUN_ROOT}/config/chunk-manifest.tsv
{
    printf 'chunk_id\tscenario_id\n'
    awk -F '\t' -v size="${SCENARIOS_PER_ARRAY_TASK}" \
        'NR > 1 {print int((NR - 2) / size) + 1 "\t" $1}' "${MANIFEST}"
} > "${CHUNK_MANIFEST}"
ARRAY_TASKS=$(( (TASKS + SCENARIOS_PER_ARRAY_TASK - 1) / SCENARIOS_PER_ARRAY_TASK ))

git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_ROOT}/provenance/git-commit.txt"
{
    printf 'field\tvalue\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'config\t%s\n' "${CONFIG}"
    printf 'acceptance_criteria\t%s\n' "${RUN_ROOT}/config/acceptance-criteria.tsv"
    printf 'n_scenarios\t%s\n' "${TASKS}"
    printf 'array_tasks\t%s\n' "${ARRAY_TASKS}"
    printf 'scenarios_per_array_task\t%s\n' "${SCENARIOS_PER_ARRAY_TASK}"
    printf 'max_concurrent\t%s\n' "${MAX_CONCURRENT}"
    printf 'gemma\t%s\n' "$(command -v gemma)"
    printf 'purpose\tbslmm_full_module02_validation_grid\n'
    printf 'observed_production_changed\tFALSE\n'
} > "${RUN_ROOT}/provenance/run-metadata.tsv"

if [[ "${DRY_RUN,,}" == "true" ]]; then
    echo "Dry run prepared ${TASKS} scenarios / ${ARRAY_TASKS} array tasks"
    echo "Run directory: ${RUN_ROOT}"
    exit 0
fi

VAL_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --partition=short \
    --array="1-${ARRAY_TASKS}%${MAX_CONCURRENT}" \
    --job-name="bslmm_val" \
    --output="${RUN_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_ENV="${ENV_PATH}",BSLMM_VAL_MANIFEST="${MANIFEST}",BSLMM_VAL_CHUNK_MANIFEST="${CHUNK_MANIFEST}",BSLMM_VAL_OUTPUT_ROOT="${RUN_ROOT}",BSLMM_VAL_WORK_ROOT="${RUN_ROOT}/work",BSLMM_VAL_KEEP_WORK="${KEEP_WORK}" \
    "${RUN_SCRIPT_DIR}/step_bslmm_validation_chunk.sh")
VAL_JOB_ID=${VAL_JOB%%;*}

EVAL_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --partition=short \
    --dependency="afterok:${VAL_JOB_ID}" \
    --job-name="bslmm_val_eval" \
    --ntasks=1 --cpus-per-task=1 --mem=24G --time=01:00:00 \
    --output="${RUN_ROOT}/logs/%x.%j.log" \
    --wrap="export PATH=/projects/p32505/opt/bin:\$PATH; ${ENV_PATH}/bin/Rscript ${RUN_SCRIPT_DIR}/17_evaluate_bslmm_validation.R --input_dir=${RUN_ROOT}/raw --manifest=${MANIFEST} --config=${RUN_ROOT}/config/bslmm-validation.tsv --criteria=${RUN_ROOT}/config/acceptance-criteria.tsv --output_dir=${RUN_ROOT}/combined --fail_on_rejection=TRUE")
EVAL_JOB_ID=${EVAL_JOB%%;*}

{
    printf 'stage\tjob_id\n'
    printf 'validation_array\t%s\n' "${VAL_JOB_ID}"
    printf 'evaluate\t%s\n' "${EVAL_JOB_ID}"
} > "${RUN_ROOT}/provenance/submitted-jobs.tsv"

echo "Submitted BSLMM Module 02 validation"
echo "  run:       ${RUN_ROOT}"
echo "  scenarios: ${TASKS}"
echo "  array:     ${VAL_JOB_ID} (${ARRAY_TASKS} tasks)"
echo "  evaluate:  ${EVAL_JOB_ID}"
echo "  decision:  ${RUN_ROOT}/combined/bslmm-validation-decision.tsv"
echo "Observed production is unchanged until decision is PASS and PI authorizes replacement."
