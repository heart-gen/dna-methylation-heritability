#!/bin/bash
#
# Paired BSLMM vs calibrated elastic-net estimator screen.
# Does not modify production Module 02 observed outputs.
#
# Usage, from the module _m/ directory:
#   ../_h/submit_bslmm_pilot.sh [RUN_ID] [CONFIG_TSV]
#
# Environment:
#   CAL_H2_ENV, SBATCH_ACCOUNT, MAX_CONCURRENT, SCENARIOS_PER_ARRAY_TASK
#   SUBMIT_BSLMM_PILOT_DRY_RUN=TRUE  prepare only
#   BSLMM_PILOT_KEEP_WORK=TRUE       retain per-scenario GEMMA work dirs

set -euo pipefail

usage() {
    echo "Usage: $0 [RUN_ID] [CONFIG_TSV]"
    echo "Default config: config/bslmm-pilot.tsv"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
if (( $# > 2 )); then
    usage >&2
    exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
DRY_RUN=${SUBMIT_BSLMM_PILOT_DRY_RUN:-FALSE}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
KEEP_WORK=${BSLMM_PILOT_KEEP_WORK:-FALSE}

CONFIG=${2:-${ANALYSIS_DIR}/config/bslmm-pilot.tsv}
CONFIG=$(readlink -f "${CONFIG}")
if [[ ! -f "${CONFIG}" ]]; then
    echo "Missing config: ${CONFIG}" >&2
    exit 1
fi

RUN_ID=${1:-bslmm-en-pilot-$(date -u +%Y%m%d)}
if [[ ! "${RUN_ID}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "Invalid RUN_ID: ${RUN_ID}" >&2
    exit 1
fi

if [[ ! -x "${ENV_PATH}/bin/Rscript" ]]; then
    echo "Missing Rscript in ${ENV_PATH}" >&2
    exit 1
fi
export PATH="/projects/p32505/opt/bin:${PATH}"
if [[ ! -x "$(command -v gemma)" ]]; then
    echo "gemma is not executable on PATH" >&2
    exit 1
fi
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

# Read defaults from config unless overridden by environment.
MAX_CONCURRENT=${MAX_CONCURRENT:-$(awk -F '\t' '$1=="max_concurrent"{print $2; exit}' "${CONFIG}")}
SCENARIOS_PER_ARRAY_TASK=${SCENARIOS_PER_ARRAY_TASK:-$(awk -F '\t' '$1=="scenarios_per_array_task"{print $2; exit}' "${CONFIG}")}
CALIBRATION_REL=$(awk -F '\t' '$1=="calibration_model_relpath"{print $2; exit}' "${CONFIG}")
EXPECTED_SHA=$(awk -F '\t' '$1=="expected_calibration_sha256"{print $2; exit}' "${CONFIG}")
CALIBRATION_MODEL=$(readlink -f "${ANALYSIS_DIR}/${CALIBRATION_REL}")
if [[ ! -f "${CALIBRATION_MODEL}" ]]; then
    echo "Missing frozen calibration model: ${CALIBRATION_MODEL}" >&2
    exit 1
fi
OBS_SHA=$(sha256sum "${CALIBRATION_MODEL}" | awk '{print $1}')
if [[ "${OBS_SHA}" != "${EXPECTED_SHA}" ]]; then
    echo "Calibration SHA mismatch: expected ${EXPECTED_SHA} got ${OBS_SHA}" >&2
    exit 1
fi

RUN_ROOT=${RUN_BASE}/${RUN_ID}
if [[ -e "${RUN_ROOT}" ]]; then
    echo "Run directory already exists: ${RUN_ROOT}" >&2
    exit 1
fi
mkdir -p "${RUN_ROOT}/config" "${RUN_ROOT}/logs" "${RUN_ROOT}/provenance" \
    "${RUN_ROOT}/code" "${RUN_ROOT}/raw" "${RUN_ROOT}/combined" "${RUN_ROOT}/work"

cp "${CONFIG}" "${RUN_ROOT}/config/bslmm-pilot.tsv"
cp -a "${SCRIPT_DIR}" "${RUN_ROOT}/code/_h"
RUN_SCRIPT_DIR=${RUN_ROOT}/code/_h
MANIFEST=${RUN_ROOT}/config/scenarios.tsv

cd "${REPO_ROOT}"
"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/12_make_bslmm_pilot_manifest.R" \
    --config="${RUN_ROOT}/config/bslmm-pilot.tsv" \
    --output="${MANIFEST}"
TASKS=$(($(wc -l < "${MANIFEST}") - 1))
if (( TASKS < 1 )); then
    echo "Manifest contains no scenarios" >&2
    exit 1
fi

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
    printf 'calibration_model\t%s\n' "${CALIBRATION_MODEL}"
    printf 'calibration_sha256\t%s\n' "${EXPECTED_SHA}"
    printf 'n_scenarios\t%s\n' "${TASKS}"
    printf 'array_tasks\t%s\n' "${ARRAY_TASKS}"
    printf 'scenarios_per_array_task\t%s\n' "${SCENARIOS_PER_ARRAY_TASK}"
    printf 'max_concurrent\t%s\n' "${MAX_CONCURRENT}"
    printf 'gemma\t%s\n' "$(command -v gemma)"
    printf 'purpose\tpaired_bslmm_vs_en_estimator_screen\n'
} > "${RUN_ROOT}/provenance/run-metadata.tsv"

if [[ "${DRY_RUN,,}" == "true" ]]; then
    echo "Dry run prepared ${TASKS} scenarios / ${ARRAY_TASKS} array tasks"
    echo "Run directory: ${RUN_ROOT}"
    exit 0
fi

PILOT_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --partition=short \
    --array="1-${ARRAY_TASKS}%${MAX_CONCURRENT}" \
    --job-name="bslmm_en_pilot" \
    --output="${RUN_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_ENV="${ENV_PATH}",BSLMM_PILOT_MANIFEST="${MANIFEST}",BSLMM_PILOT_CHUNK_MANIFEST="${CHUNK_MANIFEST}",BSLMM_PILOT_OUTPUT_ROOT="${RUN_ROOT}",BSLMM_PILOT_CALIBRATION_MODEL="${CALIBRATION_MODEL}",BSLMM_PILOT_CALIBRATION_SHA256="${EXPECTED_SHA}",BSLMM_PILOT_WORK_ROOT="${RUN_ROOT}/work",BSLMM_PILOT_KEEP_WORK="${KEEP_WORK}" \
    "${RUN_SCRIPT_DIR}/step_bslmm_pilot_chunk.sh")
PILOT_JOB_ID=${PILOT_JOB%%;*}

# Summarize via wrap so BASH_SOURCE path resolution is not broken by Slurm spool copies.
SUM_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --partition=short \
    --dependency="afterok:${PILOT_JOB_ID}" \
    --job-name="bslmm_en_sum" \
    --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00 \
    --output="${RUN_ROOT}/logs/%x.%j.log" \
    --wrap="export PATH=/projects/p32505/opt/bin:\$PATH; ${ENV_PATH}/bin/Rscript ${RUN_SCRIPT_DIR}/14_summarize_bslmm_pilot.R --input_dir=${RUN_ROOT}/raw --manifest=${MANIFEST} --config=${RUN_ROOT}/config/bslmm-pilot.tsv --output_dir=${RUN_ROOT}/combined")
SUM_JOB_ID=${SUM_JOB%%;*}

{
    printf 'stage\tjob_id\n'
    printf 'pilot_array\t%s\n' "${PILOT_JOB_ID}"
    printf 'summarize\t%s\n' "${SUM_JOB_ID}"
} > "${RUN_ROOT}/provenance/submitted-jobs.tsv"

echo "Submitted paired BSLMM-EN pilot"
echo "  run:        ${RUN_ROOT}"
echo "  scenarios:  ${TASKS}"
echo "  array job:  ${PILOT_JOB_ID} (${ARRAY_TASKS} tasks)"
echo "  summarize:  ${SUM_JOB_ID}"
echo "Decision table will be: ${RUN_ROOT}/combined/bslmm-en-pilot-decision.tsv"
