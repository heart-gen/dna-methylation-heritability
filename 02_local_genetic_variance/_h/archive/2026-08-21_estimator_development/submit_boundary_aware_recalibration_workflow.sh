#!/bin/bash

## Fresh calibration-training sims → boundary-aware/two-part candidate fit →
## freeze → entirely new validation set.
##
## The motivating failed validation grid
## (lgv-calibration-all-20260817/evaluation) is never reused.

set -euo pipefail

usage() {
    echo "Usage: $0 TRAIN_RUN_ID VALIDATION_RUN_ID [TRAIN_SEED_OFFSET] [VAL_SEED_OFFSET]"
    echo "Environment: CAL_H2_ENV, SBATCH_ACCOUNT, MAX_CONCURRENT, SCENARIOS_PER_ARRAY_TASK"
    echo "Set SUBMIT_CAL_H2_DRY_RUN=TRUE to prepare provenance without sbatch."
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
if (( $# < 2 || $# > 4 )); then
    usage >&2
    exit 1
fi

TRAIN_RUN_ID=$1
VALIDATION_RUN_ID=$2
TRAIN_SEED_OFFSET=${3:-300000000}
VAL_SEED_OFFSET=${4:-400000000}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
MAX_CONCURRENT=${MAX_CONCURRENT:-50}
SCENARIOS_PER_ARRAY_TASK=${SCENARIOS_PER_ARRAY_TASK:-10}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
DRY_RUN=${SUBMIT_CAL_H2_DRY_RUN:-FALSE}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
CONFIG=${ANALYSIS_DIR}/config/analysis.tsv
CANDIDATES=${ANALYSIS_DIR}/config/boundary-aware-candidates.tsv
PROTECTED_VALIDATION=${RUN_BASE}/lgv-calibration-all-20260817

for id in "${TRAIN_RUN_ID}" "${VALIDATION_RUN_ID}"; do
    if [[ ! "${id}" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "Run IDs may contain only letters, numbers, periods, underscores, and hyphens" >&2
        exit 1
    fi
done
for offset in "${TRAIN_SEED_OFFSET}" "${VAL_SEED_OFFSET}"; do
    if [[ ! "${offset}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Seed offsets must be positive integers" >&2
        exit 1
    fi
done
if [[ "${TRAIN_SEED_OFFSET}" -eq "${VAL_SEED_OFFSET}" ]]; then
    echo "Training and validation seed offsets must differ" >&2
    exit 1
fi
if [[ "${TRAIN_SEED_OFFSET}" -eq 0 || "${VAL_SEED_OFFSET}" -eq 0 ]]; then
    echo "Refuse seed offset 0 to avoid colliding with the motivating grid" >&2
    exit 1
fi
if [[ "${TRAIN_RUN_ID}" == "lgv-calibration-all-20260817" || \
      "${VALIDATION_RUN_ID}" == "lgv-calibration-all-20260817" ]]; then
    echo "Refuse to overwrite the motivating validation run" >&2
    exit 1
fi
if [[ ! -f "${CONFIG}" || ! -f "${CANDIDATES}" ]]; then
    echo "analysis.tsv or boundary-aware-candidates.tsv is missing" >&2
    exit 1
fi
if [[ ! -x "${ENV_PATH}/bin/Rscript" ]]; then
    echo "Dedicated calibration environment missing: ${ENV_PATH}" >&2
    exit 1
fi
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

TRAIN_ROOT=${RUN_BASE}/${TRAIN_RUN_ID}
VALIDATION_ROOT=${RUN_BASE}/${VALIDATION_RUN_ID}
if [[ -e "${TRAIN_ROOT}" ]]; then
    echo "Training run already exists: ${TRAIN_ROOT}" >&2
    exit 1
fi
if [[ -e "${VALIDATION_ROOT}" ]]; then
    echo "Validation run already exists: ${VALIDATION_ROOT}" >&2
    exit 1
fi
if [[ ! -d "${PROTECTED_VALIDATION}/evaluation" ]]; then
    echo "Warning: protected motivating validation directory not found at ${PROTECTED_VALIDATION}" >&2
fi

mkdir -p "${TRAIN_ROOT}/config" "${TRAIN_ROOT}/logs" \
    "${TRAIN_ROOT}/provenance" "${TRAIN_ROOT}/code" \
    "${TRAIN_ROOT}/calibration" "${TRAIN_ROOT}/raw/calibration"

cp "${CONFIG}" "${TRAIN_ROOT}/config/analysis.tsv"
cp "${CONFIG}" "${TRAIN_ROOT}/config/source-analysis.tsv"
cp "${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
    "${TRAIN_ROOT}/config/acceptance-criteria.tsv"
cp "${CANDIDATES}" "${TRAIN_ROOT}/config/boundary-aware-candidates.tsv"
cp "${ANALYSIS_DIR}/config/environment.yml" "${TRAIN_ROOT}/config/environment.yml"
cp -a "${SCRIPT_DIR}" "${TRAIN_ROOT}/code/_h"
RUN_SCRIPT_DIR=${TRAIN_ROOT}/code/_h

MANIFEST=${TRAIN_ROOT}/config/scenarios.tsv
"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/00_make_manifest.R" \
    --config="${TRAIN_ROOT}/config/analysis.tsv" \
    --output="${MANIFEST}" \
    --split=calibration \
    --seed-offset="${TRAIN_SEED_OFFSET}"
TASKS=$(($(wc -l < "${MANIFEST}") - 1))
if (( TASKS < 1 )); then
    echo "Training manifest contains no tasks" >&2
    exit 1
fi
SIMULATION_CHUNK_MANIFEST=${TRAIN_ROOT}/config/simulation-chunk-manifest.tsv
{
    printf 'chunk_id\tscenario_id\n'
    awk -F '\t' -v size="${SCENARIOS_PER_ARRAY_TASK}" \
        'NR > 1 {print int((NR - 2) / size) + 1 "\t" $1}' "${MANIFEST}"
} > "${SIMULATION_CHUNK_MANIFEST}"
SIMULATION_JOBS=$(( (TASKS + SCENARIOS_PER_ARRAY_TASK - 1) / SCENARIOS_PER_ARRAY_TASK ))

{
    printf 'field\tvalue\n'
    printf 'train_run_id\t%s\n' "${TRAIN_RUN_ID}"
    printf 'validation_run_id\t%s\n' "${VALIDATION_RUN_ID}"
    printf 'train_seed_offset\t%s\n' "${TRAIN_SEED_OFFSET}"
    printf 'validation_seed_offset\t%s\n' "${VAL_SEED_OFFSET}"
    printf 'calibration_tasks\t%s\n' "${TASKS}"
    printf 'simulation_array_tasks\t%s\n' "${SIMULATION_JOBS}"
    printf 'protected_validation_run\tlgv-calibration-all-20260817\n'
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "${TRAIN_ROOT}/provenance/run-metadata.tsv"
git -C "${REPO_ROOT}" rev-parse HEAD > "${TRAIN_ROOT}/provenance/git-commit.txt"
conda list -p "${ENV_PATH}" --explicit > \
    "${TRAIN_ROOT}/provenance/conda-explicit-spec.txt"
for path in \
    "${RUN_SCRIPT_DIR}"/*.R \
    "${RUN_SCRIPT_DIR}"/*.sh \
    "${TRAIN_ROOT}/config"/*; do
    sha256sum "${path}"
done > "${TRAIN_ROOT}/provenance/sha256sums.txt"

if [[ "${DRY_RUN,,}" == "true" ]]; then
    echo "Dry run prepared ${TASKS} calibration scenarios without submitting jobs"
    echo "Training run: ${TRAIN_ROOT}"
    echo "Validation run (not yet created): ${VALIDATION_ROOT}"
    exit 0
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch is unavailable; set SUBMIT_CAL_H2_DRY_RUN=TRUE to prepare only" >&2
    exit 1
fi

SIM_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --array="1-${SIMULATION_JOBS}%${MAX_CONCURRENT}" \
    --job-name="cal_h2_train" \
    --output="${TRAIN_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",SCENARIO_MANIFEST="${MANIFEST}",CAL_H2_RECOVERY_MANIFEST="${SIMULATION_CHUNK_MANIFEST}",SIMULATION_OUTPUT_ROOT="${TRAIN_ROOT}/raw" \
    "${RUN_SCRIPT_DIR}/step_2_recover_simulation.sh")
SIM_JOB_ID=${SIM_JOB%%;*}

FIT_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${SIM_JOB_ID}" \
    --job-name="cal_h2_fit_boundary" \
    --output="${TRAIN_ROOT}/logs/%x.%j.log" \
    --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_RUN_ROOT="${TRAIN_ROOT}",CAL_H2_CANDIDATES="${TRAIN_ROOT}/config/boundary-aware-candidates.tsv" \
    "${RUN_SCRIPT_DIR}/step_3b_fit_boundary_aware_calibration.sh")
FIT_JOB_ID=${FIT_JOB%%;*}

## After freeze, create an independent validation run and evaluate the frozen model.
VALIDATE_WRAPPER=${TRAIN_ROOT}/provenance/submit_fresh_validation.sh
cat > "${VALIDATE_WRAPPER}" <<EOF
#!/bin/bash
set -euo pipefail
export SBATCH_ACCOUNT=${ACCOUNT}
export MAX_CONCURRENT=${MAX_CONCURRENT}
export CAL_H2_ENV=${ENV_PATH}
export CAL_H2_RUN_BASE=${RUN_BASE}
bash "${RUN_SCRIPT_DIR}/submit_fresh_validation_workflow.sh" \
    "${TRAIN_RUN_ID}" "${VALIDATION_RUN_ID}" "${VAL_SEED_OFFSET}"
EOF
chmod +x "${VALIDATE_WRAPPER}"

VAL_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${FIT_JOB_ID}" \
    --job-name="cal_h2_launch_val" \
    --partition=short \
    --time=00:20:00 \
    --mem=2G \
    --output="${TRAIN_ROOT}/logs/%x.%j.log" \
    --wrap="bash ${VALIDATE_WRAPPER}")
VAL_JOB_ID=${VAL_JOB%%;*}

{
    printf 'stage\tjob_id\tdependency\n'
    printf 'calibration_training_simulations\t%s\tNA\n' "${SIM_JOB_ID}"
    printf 'boundary_aware_candidate_fit\t%s\tafterok:%s\n' "${FIT_JOB_ID}" "${SIM_JOB_ID}"
    printf 'launch_fresh_validation\t%s\tafterok:%s\n' "${VAL_JOB_ID}" "${FIT_JOB_ID}"
} > "${TRAIN_ROOT}/provenance/submitted-jobs.tsv"

echo "Training run: ${TRAIN_ROOT}"
echo "Calibration-only scenarios: ${TASKS} (${SIMULATION_JOBS} array tasks)"
echo "Simulation array: ${SIM_JOB_ID}"
echo "Candidate fit / freeze: ${FIT_JOB_ID}"
echo "Fresh validation launcher: ${VAL_JOB_ID} -> ${VALIDATION_RUN_ID}"
echo "Protected untouched evidence: ${PROTECTED_VALIDATION}"
