#!/bin/bash

set -euo pipefail

usage() {
    echo "Usage: $0 [RUN_ID] [CONFIG_TSV]"
    echo "Environment: CAL_H2_ENV, SBATCH_ACCOUNT, MAX_CONCURRENT"
    echo "Set SUBMIT_CAL_H2_DRY_RUN=TRUE to prepare provenance without sbatch."
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
MAX_CONCURRENT=${MAX_CONCURRENT:-200}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
RUN_ID=${1:-$(date -u +%Y%m%dT%H%M%SZ)}
CONFIG=${2:-${ANALYSIS_DIR}/config/analysis.tsv}
DRY_RUN=${SUBMIT_CAL_H2_DRY_RUN:-FALSE}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}

if [[ ! "${RUN_ID}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "RUN_ID may contain only letters, numbers, periods, underscores, and hyphens" >&2
    exit 1
fi
if [[ ! "${MAX_CONCURRENT}" =~ ^[1-9][0-9]*$ ]]; then
    echo "MAX_CONCURRENT must be a positive integer" >&2
    exit 1
fi
if [[ ! -f "${CONFIG}" ]]; then
    echo "Configuration file is missing: ${CONFIG}" >&2
    exit 1
fi
CONFIG=$(readlink -f "${CONFIG}")
if [[ ! -x "${ENV_PATH}/bin/Rscript" ]]; then
    echo "The dedicated calibration environment is missing or incomplete: ${ENV_PATH}" >&2
    echo "Create it with: bash ${SCRIPT_DIR}/setup_environment.sh" >&2
    exit 1
fi
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

RUN_ROOT=${RUN_BASE}/${RUN_ID}
if [[ "${RUN_ROOT}" == *","* ]]; then
    echo "Run paths containing commas cannot be exported safely through sbatch" >&2
    exit 1
fi
if [[ -e "${RUN_ROOT}" ]]; then
    echo "Run directory already exists; choose a new RUN_ID: ${RUN_ROOT}" >&2
    exit 1
fi
mkdir -p "${RUN_ROOT}/config" "${RUN_ROOT}/logs" \
    "${RUN_ROOT}/provenance" "${RUN_ROOT}/code"

cp "${CONFIG}" "${RUN_ROOT}/config/analysis.tsv"
cp "${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
    "${RUN_ROOT}/config/acceptance-criteria.tsv"
cp "${ANALYSIS_DIR}/config/environment.yml" "${RUN_ROOT}/config/environment.yml"
cp -a "${SCRIPT_DIR}" "${RUN_ROOT}/code/_h"
RUN_SCRIPT_DIR=${RUN_ROOT}/code/_h
MANIFEST=${RUN_ROOT}/config/scenarios.tsv

cd "${REPO_ROOT}"
"${RUN_SCRIPT_DIR}/step_1_generate_manifest.sh" \
    "${RUN_ROOT}/config/analysis.tsv" "${MANIFEST}"
TASKS=$(($(wc -l < "${MANIFEST}") - 1))
if (( TASKS < 1 )); then
    echo "Manifest contains no tasks" >&2
    exit 1
fi

git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_ROOT}/provenance/git-commit.txt"
conda list -p "${ENV_PATH}" --explicit > \
    "${RUN_ROOT}/provenance/conda-explicit-spec.txt"
{
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'repo_root\t%s\n' "${REPO_ROOT}"
    printf 'source_config\t%s\n' "${CONFIG}"
    printf 'environment\t%s\n' "${ENV_PATH}"
    printf 'account\t%s\n' "${ACCOUNT}"
    printf 'max_concurrent\t%s\n' "${MAX_CONCURRENT}"
    printf 'scenario_tasks\t%s\n' "${TASKS}"
} > "${RUN_ROOT}/provenance/run-metadata.tsv"
for path in \
    "${RUN_SCRIPT_DIR}"/*.R \
    "${RUN_SCRIPT_DIR}"/*.sh \
    "${RUN_ROOT}/config"/*; do
    sha256sum "${path}"
done > "${RUN_ROOT}/provenance/sha256sums.txt"

if [[ "${DRY_RUN,,}" == "true" ]]; then
    echo "Dry run prepared ${TASKS} scenarios without submitting jobs"
    echo "Run directory: ${RUN_ROOT}"
    exit 0
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch is unavailable; set SUBMIT_CAL_H2_DRY_RUN=TRUE to prepare only" >&2
    exit 1
fi

SIM_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --array="1-${TASKS}%${MAX_CONCURRENT}" \
    --output="${RUN_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_RUN_ROOT="${RUN_ROOT}",SCENARIO_MANIFEST="${MANIFEST}",SIMULATION_OUTPUT_ROOT="${RUN_ROOT}/raw" \
    "${RUN_SCRIPT_DIR}/step_2_simulate_and_crossfit.sh")
SIM_JOB_ID=${SIM_JOB%%;*}
CAL_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${SIM_JOB_ID}" \
    --output="${RUN_ROOT}/logs/%x.%j.log" \
    --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_RUN_ROOT="${RUN_ROOT}" \
    "${RUN_SCRIPT_DIR}/step_3_fit_calibration.sh")
CAL_JOB_ID=${CAL_JOB%%;*}
EVAL_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${CAL_JOB_ID}" \
    --output="${RUN_ROOT}/logs/%x.%j.log" \
    --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_RUN_ROOT="${RUN_ROOT}",CAL_H2_FAIL_ON_REJECTION=TRUE \
    "${RUN_SCRIPT_DIR}/step_4_evaluate_calibration.sh")
EVAL_JOB_ID=${EVAL_JOB%%;*}

{
    printf 'stage\tjob_id\tdependency\n'
    printf 'simulation\t%s\tNA\n' "${SIM_JOB_ID}"
    printf 'calibration\t%s\tafterok:%s\n' "${CAL_JOB_ID}" "${SIM_JOB_ID}"
    printf 'evaluation\t%s\tafterok:%s\n' "${EVAL_JOB_ID}" "${CAL_JOB_ID}"
} > "${RUN_ROOT}/provenance/submitted-jobs.tsv"

echo "Run directory: ${RUN_ROOT}"
echo "Simulation array: ${SIM_JOB_ID}"
echo "Calibration: ${CAL_JOB_ID}"
echo "Evaluation and acceptance gate: ${EVAL_JOB_ID}"
