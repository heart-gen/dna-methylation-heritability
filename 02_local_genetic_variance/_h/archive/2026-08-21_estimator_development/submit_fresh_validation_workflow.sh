#!/bin/bash

set -euo pipefail

usage() {
    echo "Usage: $0 MODEL_RUN_ID VALIDATION_RUN_ID [SEED_OFFSET]"
    echo "Runs a new evaluation split against a frozen calibration model."
}
if (( $# < 2 || $# > 3 )); then
    usage >&2
    exit 1
fi

MODEL_RUN_ID=$1
VALIDATION_RUN_ID=$2
SEED_OFFSET=${3:-100000000}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
MAX_CONCURRENT=${MAX_CONCURRENT:-50}
SCENARIOS_PER_ARRAY_TASK=${SCENARIOS_PER_ARRAY_TASK:-10}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
MODEL_ROOT=${RUN_BASE}/${MODEL_RUN_ID}
VALIDATION_ROOT=${RUN_BASE}/${VALIDATION_RUN_ID}
MODEL=${MODEL_ROOT}/calibration/elastic-net-calibration.rds
SOURCE_CONFIG=${MODEL_ROOT}/config/source-analysis.tsv
if [[ ! -s "${SOURCE_CONFIG}" ]]; then
    SOURCE_CONFIG=${MODEL_ROOT}/config/analysis.tsv
fi

for id in "${MODEL_RUN_ID}" "${VALIDATION_RUN_ID}"; do
    if [[ ! "${id}" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "Run IDs contain unsupported characters: ${id}" >&2
        exit 1
    fi
done
if [[ ! "${SEED_OFFSET}" =~ ^[1-9][0-9]*$ ]]; then
    echo "SEED_OFFSET must be a positive integer" >&2
    exit 1
fi
if [[ "${VALIDATION_RUN_ID}" == "lgv-calibration-all-20260817" ]]; then
    echo "Refuse to overwrite the motivating validation run" >&2
    exit 1
fi
if [[ ! -s "${MODEL}" || ! -s "${SOURCE_CONFIG}" ]]; then
    echo "Frozen model or its source configuration is missing" >&2
    exit 1
fi
if [[ -e "${VALIDATION_ROOT}" ]]; then
    echo "Validation run already exists: ${VALIDATION_ROOT}" >&2
    exit 1
fi
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

mkdir -p "${VALIDATION_ROOT}/config" "${VALIDATION_ROOT}/code" \
    "${VALIDATION_ROOT}/calibration" "${VALIDATION_ROOT}/logs" \
    "${VALIDATION_ROOT}/provenance" "${VALIDATION_ROOT}/raw/evaluation"
cp -a "${SCRIPT_DIR}" "${VALIDATION_ROOT}/code/_h"
RUN_SCRIPT_DIR=${VALIDATION_ROOT}/code/_h
cp "${SOURCE_CONFIG}" "${VALIDATION_ROOT}/config/analysis.tsv"
cp "${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
    "${VALIDATION_ROOT}/config/acceptance-criteria.tsv"
if [[ -f "${ANALYSIS_DIR}/config/environment.yml" ]]; then
    cp "${ANALYSIS_DIR}/config/environment.yml" "${VALIDATION_ROOT}/config/environment.yml"
elif [[ -f "${ANALYSIS_DIR}/environment.yml" ]]; then
    cp "${ANALYSIS_DIR}/environment.yml" "${VALIDATION_ROOT}/config/environment.yml"
fi
cp "${MODEL}" "${VALIDATION_ROOT}/calibration/elastic-net-calibration.rds"
SELECTED="${MODEL_ROOT}/calibration/selected-candidate.tsv"
if [[ -s "${SELECTED}" ]]; then
    cp "${SELECTED}" "${VALIDATION_ROOT}/calibration/selected-candidate.tsv"
fi

MANIFEST=${VALIDATION_ROOT}/config/fresh-evaluation-scenarios.tsv
"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/00_make_manifest.R" \
    --config="${VALIDATION_ROOT}/config/analysis.tsv" \
    --output="${MANIFEST}" \
    --split=evaluation \
    --seed-offset="${SEED_OFFSET}"
TASKS=$(($(wc -l < "${MANIFEST}") - 1))
if (( TASKS < 1 )); then
    echo "Fresh evaluation manifest contains no tasks" >&2
    exit 1
fi
SIMULATION_CHUNK_MANIFEST=${VALIDATION_ROOT}/config/simulation-chunk-manifest.tsv
{
    printf 'chunk_id\tscenario_id\n'
    awk -F '\t' -v size="${SCENARIOS_PER_ARRAY_TASK}" \
        'NR > 1 {print int((NR - 2) / size) + 1 "\t" $1}' "${MANIFEST}"
} > "${SIMULATION_CHUNK_MANIFEST}"
SIMULATION_JOBS=$(( (TASKS + SCENARIOS_PER_ARRAY_TASK - 1) / SCENARIOS_PER_ARRAY_TASK ))

{
    printf 'field\tvalue\n'
    printf 'model_run_id\t%s\n' "${MODEL_RUN_ID}"
    printf 'validation_run_id\t%s\n' "${VALIDATION_RUN_ID}"
    printf 'seed_offset\t%s\n' "${SEED_OFFSET}"
    printf 'evaluation_tasks\t%s\n' "${TASKS}"
    printf 'simulation_array_tasks\t%s\n' "${SIMULATION_JOBS}"
    printf 'scenarios_per_array_task\t%s\n' "${SCENARIOS_PER_ARRAY_TASK}"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'model_sha256\t%s\n' "$(sha256sum "${MODEL}" | awk '{print $1}')"
    printf 'protected_validation_run\tlgv-calibration-all-20260817\n'
} > "${VALIDATION_ROOT}/provenance/validation-metadata.tsv"
git -C "${REPO_ROOT}" rev-parse HEAD > "${VALIDATION_ROOT}/provenance/git-commit.txt"
conda list -p "${ENV_PATH}" --explicit > \
    "${VALIDATION_ROOT}/provenance/conda-explicit-spec.txt"
for path in "${RUN_SCRIPT_DIR}"/*.R "${RUN_SCRIPT_DIR}"/*.sh \
    "${VALIDATION_ROOT}/config"/* "${VALIDATION_ROOT}/calibration"/*; do
    sha256sum "${path}"
done > "${VALIDATION_ROOT}/provenance/sha256sums.txt"

cd "${REPO_ROOT}"
SIM_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --array="1-${SIMULATION_JOBS}%${MAX_CONCURRENT}" \
    --job-name="cal_h2_val" \
    --output="${VALIDATION_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",SCENARIO_MANIFEST="${MANIFEST}",CAL_H2_RECOVERY_MANIFEST="${SIMULATION_CHUNK_MANIFEST}",SIMULATION_OUTPUT_ROOT="${VALIDATION_ROOT}/raw" \
    "${RUN_SCRIPT_DIR}/step_2_recover_simulation.sh")
SIM_JOB_ID=${SIM_JOB%%;*}
EVAL_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${SIM_JOB_ID}" \
    --job-name="cal_h2_eval" \
    --output="${VALIDATION_ROOT}/logs/%x.%j.log" \
    --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_RUN_ROOT="${VALIDATION_ROOT}",CAL_H2_FAIL_ON_REJECTION=TRUE \
    "${RUN_SCRIPT_DIR}/step_4_evaluate_calibration.sh")
EVAL_JOB_ID=${EVAL_JOB%%;*}
{
    printf 'stage\tjob_id\tdependency\n'
    printf 'fresh_evaluation_simulations\t%s\tNA\n' "${SIM_JOB_ID}"
    printf 'frozen_model_evaluation\t%s\tafterok:%s\n' "${EVAL_JOB_ID}" "${SIM_JOB_ID}"
} > "${VALIDATION_ROOT}/provenance/submitted-jobs.tsv"

echo "Validation run: ${VALIDATION_ROOT}"
echo "Fresh evaluation array: ${SIM_JOB_ID} (${SIMULATION_JOBS} jobs for ${TASKS} scenarios)"
echo "Frozen-model evaluation gate: ${EVAL_JOB_ID}"
