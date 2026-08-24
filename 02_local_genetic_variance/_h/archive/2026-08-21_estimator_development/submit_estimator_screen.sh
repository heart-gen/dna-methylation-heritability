#!/bin/bash

## Paired estimator-settings screen: do more outer repeats, lambda.min, or a
## ridge-ward alpha grid reduce raw-statistic dispersion?
##
## Screening only. Writes no calibration model, changes no production output,
## and does not recalibrate Module 02 (AGENTS.md 7.2).

set -euo pipefail

usage() {
    echo "Usage: $0 RUN_ID [CONFIG_TSV] [ARMS_TSV]"
    echo "Environment: CAL_H2_ENV, SBATCH_ACCOUNT, MAX_CONCURRENT,"
    echo "             SCENARIOS_PER_ARRAY_TASK, SCREEN_PARTITION, SCREEN_TIME"
    echo "Set SUBMIT_CAL_H2_DRY_RUN=TRUE to prepare provenance without sbatch."
}
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
if (( $# < 1 || $# > 3 )); then
    usage >&2
    exit 1
fi

RUN_ID=$1
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
CONFIG=${2:-${ANALYSIS_DIR}/config/estimator-screen-20260818.tsv}
ARMS=${3:-${ANALYSIS_DIR}/config/estimator-screen-arms.tsv}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
MAX_CONCURRENT=${MAX_CONCURRENT:-50}
SCENARIOS_PER_ARRAY_TASK=${SCENARIOS_PER_ARRAY_TASK:-4}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
DRY_RUN=${SUBMIT_CAL_H2_DRY_RUN:-FALSE}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
## Heavy arms triple the repeats and double the alpha grid, and the runtime
## distribution has a long tail, so default to a partition that allows >4 h.
SCREEN_PARTITION=${SCREEN_PARTITION:-normal}
SCREEN_TIME=${SCREEN_TIME:-16:00:00}

if [[ ! "${RUN_ID}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "RUN_ID may contain only letters, numbers, periods, underscores, and hyphens" >&2
    exit 1
fi
if [[ ! -f "${CONFIG}" || ! -f "${ARMS}" ]]; then
    echo "Screen configuration or arm table is missing" >&2
    exit 1
fi
if [[ ! -x "${ENV_PATH}/bin/Rscript" ]]; then
    echo "The dedicated calibration environment is missing: ${ENV_PATH}" >&2
    exit 1
fi
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

RUN_ROOT=${RUN_BASE}/${RUN_ID}
if [[ -e "${RUN_ROOT}" ]]; then
    echo "Run directory already exists; choose a new RUN_ID: ${RUN_ROOT}" >&2
    exit 1
fi
mkdir -p "${RUN_ROOT}/config" "${RUN_ROOT}/logs" "${RUN_ROOT}/provenance" \
    "${RUN_ROOT}/code" "${RUN_ROOT}/raw" "${RUN_ROOT}/summary"
cp "${CONFIG}" "${RUN_ROOT}/config/estimator-screen.tsv"
cp "${ARMS}" "${RUN_ROOT}/config/estimator-screen-arms.tsv"
cp -a "${SCRIPT_DIR}" "${RUN_ROOT}/code/_h"
RUN_SCRIPT_DIR=${RUN_ROOT}/code/_h

MANIFEST=${RUN_ROOT}/config/screen-scenarios.tsv
"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/19_make_estimator_screen_manifest.R" \
    --config="${RUN_ROOT}/config/estimator-screen.tsv" \
    --arms="${RUN_ROOT}/config/estimator-screen-arms.tsv" \
    --output="${MANIFEST}"
TASKS=$(($(wc -l < "${MANIFEST}") - 1))
if (( TASKS < 1 )); then
    echo "Screen manifest contains no tasks" >&2
    exit 1
fi
CHUNK_MANIFEST=${RUN_ROOT}/config/screen-chunk-manifest.tsv
{
    printf 'chunk_id\tscenario_id\n'
    awk -F '\t' -v size="${SCENARIOS_PER_ARRAY_TASK}" \
        'NR > 1 {print int((NR - 2) / size) + 1 "\t" $1}' "${MANIFEST}"
} > "${CHUNK_MANIFEST}"
CHUNKS=$(( (TASKS + SCENARIOS_PER_ARRAY_TASK - 1) / SCENARIOS_PER_ARRAY_TASK ))

git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_ROOT}/provenance/git-commit.txt"
conda list -p "${ENV_PATH}" --explicit > \
    "${RUN_ROOT}/provenance/conda-explicit-spec.txt"
{
    printf 'field\tvalue\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'purpose\testimator_screening_only\n'
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'scenarios\t%s\n' "${TASKS}"
    printf 'array_tasks\t%s\n' "${CHUNKS}"
    printf 'scenarios_per_array_task\t%s\n' "${SCENARIOS_PER_ARRAY_TASK}"
    printf 'partition\t%s\n' "${SCREEN_PARTITION}"
    printf 'time_limit\t%s\n' "${SCREEN_TIME}"
    printf 'changes_production\tFALSE\n'
} > "${RUN_ROOT}/provenance/run-metadata.tsv"
for path in "${RUN_SCRIPT_DIR}"/*.R "${RUN_SCRIPT_DIR}"/*.sh \
    "${RUN_ROOT}/config"/*; do
    sha256sum "${path}"
done > "${RUN_ROOT}/provenance/sha256sums.txt"

if [[ "${DRY_RUN,,}" == "true" ]]; then
    echo "Dry run prepared ${TASKS} scenarios in ${CHUNKS} array tasks"
    echo "Run directory: ${RUN_ROOT}"
    exit 0
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch is unavailable; set SUBMIT_CAL_H2_DRY_RUN=TRUE to prepare only" >&2
    exit 1
fi

cd "${REPO_ROOT}"
SIM_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --partition="${SCREEN_PARTITION}" \
    --time="${SCREEN_TIME}" \
    --array="1-${CHUNKS}%${MAX_CONCURRENT}" \
    --job-name="cal_h2_screen" \
    --output="${RUN_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",SCENARIO_MANIFEST="${MANIFEST}",CAL_H2_RECOVERY_MANIFEST="${CHUNK_MANIFEST}",SIMULATION_OUTPUT_ROOT="${RUN_ROOT}/raw",CAL_H2_ENV="${ENV_PATH}" \
    "${RUN_SCRIPT_DIR}/step_estimator_screen_chunk.sh")
SIM_JOB_ID=${SIM_JOB%%;*}

SUM_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${SIM_JOB_ID}" \
    --partition=short \
    --time=01:00:00 \
    --mem=16G \
    --job-name="cal_h2_screen_sum" \
    --output="${RUN_ROOT}/logs/%x.%j.log" \
    --wrap="${ENV_PATH}/bin/Rscript ${RUN_SCRIPT_DIR}/20_summarize_estimator_screen.R --input=${RUN_ROOT}/raw --manifest=${MANIFEST} --config=${RUN_ROOT}/config/estimator-screen.tsv --output-dir=${RUN_ROOT}/summary")
SUM_JOB_ID=${SUM_JOB%%;*}

{
    printf 'stage\tjob_id\tdependency\n'
    printf 'screen_simulations\t%s\tNA\n' "${SIM_JOB_ID}"
    printf 'screen_summary\t%s\tafterok:%s\n' "${SUM_JOB_ID}" "${SIM_JOB_ID}"
} > "${RUN_ROOT}/provenance/submitted-jobs.tsv"

echo "Run directory: ${RUN_ROOT}"
echo "Scenarios: ${TASKS} across ${CHUNKS} array tasks (${SCREEN_PARTITION}, ${SCREEN_TIME})"
echo "Screen array: ${SIM_JOB_ID}"
echo "Screen summary: ${SUM_JOB_ID}"
