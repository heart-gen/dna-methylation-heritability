#!/bin/bash
#
# Final Module 02 absolute-PVE experiment. The family and terminal rule are
# locked in config/FINAL_JOINT_PVE_STRATEGY.md.
#
# Usage:
#   ./submit_final_joint_pve_workflow.sh [TRAIN_RUN_ID] [VALIDATION_RUN_ID] [CONFIG]

set -euo pipefail

usage() {
    echo "Usage: $0 [TRAIN_RUN_ID] [VALIDATION_RUN_ID] [CONFIG_TSV]"
}
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# > 3 )); then usage >&2; exit 1; fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
DRY_RUN=${SUBMIT_JOINT_PVE_DRY_RUN:-FALSE}
KEEP_WORK=${JOINT_PVE_KEEP_WORK:-FALSE}
CONFIG=${3:-${ANALYSIS_DIR}/config/joint-pve-20260820.tsv}
CONFIG=$(readlink -f "${CONFIG}")
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
CRITERIA=${ANALYSIS_DIR}/config/joint-pve-acceptance-criteria.tsv
STRATEGY=${ANALYSIS_DIR}/config/FINAL_JOINT_PVE_STRATEGY.md
[[ -f "${CRITERIA}" && -f "${STRATEGY}" ]] || {
    echo "Missing final joint-PVE lock files" >&2; exit 1;
}

config_value() {
    awk -F '\t' -v key="$1" '$1==key {print $2; exit}' "${CONFIG}"
}
TRAIN_ID=${1:-$(config_value training_run_id)}
VALIDATION_ID=${2:-$(config_value validation_run_id)}
for run_id in "${TRAIN_ID}" "${VALIDATION_ID}"; do
    [[ "${run_id}" =~ ^[A-Za-z0-9._-]+$ ]] || {
        echo "Invalid run ID: ${run_id}" >&2; exit 1;
    }
done
[[ "${TRAIN_ID}" != "${VALIDATION_ID}" ]] || {
    echo "Training and validation run IDs must differ" >&2; exit 1;
}
[[ -x "${ENV_PATH}/bin/Rscript" ]] || {
    echo "Missing Rscript: ${ENV_PATH}/bin/Rscript" >&2; exit 1;
}
export PATH="/projects/p32505/opt/bin:${PATH}"
command -v gemma >/dev/null || { echo "gemma not on PATH" >&2; exit 1; }
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

DEVELOPMENT_REL=$(config_value development_bslmm_run_relpath)
DEVELOPMENT_ROOT=${ANALYSIS_DIR}/${DEVELOPMENT_REL}
DEVELOPMENT_MANIFEST=${DEVELOPMENT_ROOT}/config/scenarios.tsv
DEVELOPMENT_RAW=${DEVELOPMENT_ROOT}/raw
[[ -f "${DEVELOPMENT_MANIFEST}" && -d "${DEVELOPMENT_RAW}" ]] || {
    echo "Missing locked development BSLMM run: ${DEVELOPMENT_ROOT}" >&2; exit 1;
}
EXPECTED_DEVELOPMENT=$(($(wc -l < "${DEVELOPMENT_MANIFEST}") - 1))
COMPLETED_DEVELOPMENT=$(find "${DEVELOPMENT_RAW}" -maxdepth 1 \
    -type f -name 'scenario-*.tsv' | wc -l)
[[ "${COMPLETED_DEVELOPMENT}" -eq "${EXPECTED_DEVELOPMENT}" ]] || {
    echo "Development BSLMM run is incomplete" >&2; exit 1;
}

TRAIN_ROOT=${RUN_BASE}/${TRAIN_ID}
VALIDATION_ROOT=${RUN_BASE}/${VALIDATION_ID}
for root in "${TRAIN_ROOT}" "${VALIDATION_ROOT}"; do
    [[ ! -e "${root}" ]] || { echo "Run already exists: ${root}" >&2; exit 1; }
    mkdir -p "${root}/config" "${root}/logs" "${root}/provenance" \
        "${root}/code" "${root}/raw" "${root}/combined" "${root}/work"
    cp "${CONFIG}" "${root}/config/joint-pve-20260820.tsv"
    cp "${CRITERIA}" "${root}/config/joint-pve-acceptance-criteria.tsv"
    cp "${STRATEGY}" "${root}/config/FINAL_JOINT_PVE_STRATEGY.md"
    cp -a "${SCRIPT_DIR}" "${root}/code/_h"
done
TRAIN_SCRIPT_DIR=${TRAIN_ROOT}/code/_h
VALIDATION_SCRIPT_DIR=${VALIDATION_ROOT}/code/_h
TRAIN_MANIFEST=${TRAIN_ROOT}/config/scenarios.tsv
VALIDATION_MANIFEST=${VALIDATION_ROOT}/config/scenarios.tsv

"${ENV_PATH}/bin/Rscript" "${TRAIN_SCRIPT_DIR}/21_make_joint_pve_manifest.R" \
    --config="${TRAIN_ROOT}/config/joint-pve-20260820.tsv" \
    --mode=development \
    --development_manifest="${DEVELOPMENT_MANIFEST}" \
    --output="${TRAIN_MANIFEST}"
"${ENV_PATH}/bin/Rscript" "${VALIDATION_SCRIPT_DIR}/21_make_joint_pve_manifest.R" \
    --config="${VALIDATION_ROOT}/config/joint-pve-20260820.tsv" \
    --mode=validation \
    --output="${VALIDATION_MANIFEST}"

SCENARIOS_PER_TASK=${SCENARIOS_PER_ARRAY_TASK:-$(config_value scenarios_per_array_task)}
MAX_CONCURRENT=${MAX_CONCURRENT:-$(config_value max_concurrent)}
make_chunks() {
    local manifest=$1
    local output=$2
    {
        printf 'chunk_id\tscenario_id\n'
        awk -F '\t' -v size="${SCENARIOS_PER_TASK}" \
            'NR > 1 {print int((NR - 2) / size) + 1 "\t" $1}' "${manifest}"
    } > "${output}"
}
TRAIN_CHUNKS=${TRAIN_ROOT}/config/chunk-manifest.tsv
VALIDATION_CHUNKS=${VALIDATION_ROOT}/config/chunk-manifest.tsv
make_chunks "${TRAIN_MANIFEST}" "${TRAIN_CHUNKS}"
make_chunks "${VALIDATION_MANIFEST}" "${VALIDATION_CHUNKS}"
TRAIN_SCENARIOS=$(($(wc -l < "${TRAIN_MANIFEST}") - 1))
VALIDATION_SCENARIOS=$(($(wc -l < "${VALIDATION_MANIFEST}") - 1))
TRAIN_TASKS=$(( (TRAIN_SCENARIOS + SCENARIOS_PER_TASK - 1) / SCENARIOS_PER_TASK ))
VALIDATION_TASKS=$(( (VALIDATION_SCENARIOS + SCENARIOS_PER_TASK - 1) / SCENARIOS_PER_TASK ))

COMMIT=$(git -C "${REPO_ROOT}" rev-parse HEAD)
CONFIG_SHA=$(sha256sum "${CONFIG}" | awk '{print $1}')
CRITERIA_SHA=$(sha256sum "${CRITERIA}" | awk '{print $1}')
SOURCE_SHA=$(sha256sum "${DEVELOPMENT_MANIFEST}" | awk '{print $1}')
for root in "${TRAIN_ROOT}" "${VALIDATION_ROOT}"; do
    printf '%s\n' "${COMMIT}" > "${root}/provenance/git-commit.txt"
    {
        printf 'field\tvalue\n'
        printf 'decision_lock\t2026-08-20\n'
        printf 'config_sha256\t%s\n' "${CONFIG_SHA}"
        printf 'criteria_sha256\t%s\n' "${CRITERIA_SHA}"
        printf 'development_manifest_sha256\t%s\n' "${SOURCE_SHA}"
        printf 'development_source\t%s\n' "${DEVELOPMENT_ROOT}"
        printf 'gemma\t%s\n' "$(command -v gemma)"
        printf 'terminal_experiment\tTRUE\n'
    } > "${root}/provenance/run-metadata.tsv"
done
{
    printf 'run_id\t%s\n' "${TRAIN_ID}"
    printf 'purpose\tjoint_pve_development_feature_augmentation_and_fit\n'
    printf 'n_scenarios\t%s\n' "${TRAIN_SCENARIOS}"
} >> "${TRAIN_ROOT}/provenance/run-metadata.tsv"
{
    printf 'run_id\t%s\n' "${VALIDATION_ID}"
    printf 'purpose\tuntouched_independent_joint_pve_acceptance\n'
    printf 'n_scenarios\t%s\n' "${VALIDATION_SCENARIOS}"
} >> "${VALIDATION_ROOT}/provenance/run-metadata.tsv"

if [[ "${DRY_RUN,,}" == "true" ]]; then
    echo "Dry run prepared final joint-PVE experiment"
    echo "  development: ${TRAIN_SCENARIOS} scenarios / ${TRAIN_TASKS} tasks"
    echo "  validation:  ${VALIDATION_SCENARIOS} scenarios / ${VALIDATION_TASKS} tasks"
    echo "  train root:  ${TRAIN_ROOT}"
    echo "  validation:  ${VALIDATION_ROOT}"
    exit 0
fi

TRAIN_JOB=$(sbatch --parsable --account="${ACCOUNT}" --partition=short \
    --array="1-${TRAIN_TASKS}%${MAX_CONCURRENT}" --job-name=joint_pve_train \
    --output="${TRAIN_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,JOINT_PVE_SCRIPT_DIR="${TRAIN_SCRIPT_DIR}",CAL_H2_ENV="${ENV_PATH}",JOINT_PVE_MANIFEST="${TRAIN_MANIFEST}",JOINT_PVE_CHUNK_MANIFEST="${TRAIN_CHUNKS}",JOINT_PVE_OUTPUT_ROOT="${TRAIN_ROOT}",JOINT_PVE_WORK_ROOT="${TRAIN_ROOT}/work",JOINT_PVE_DEVELOPMENT_BSLMM_ROOT="${DEVELOPMENT_RAW}",JOINT_PVE_KEEP_WORK="${KEEP_WORK}" \
    "${TRAIN_SCRIPT_DIR}/step_joint_pve_chunk.sh")
TRAIN_JOB_ID=${TRAIN_JOB%%;*}

VALIDATION_JOB=$(sbatch --parsable --account="${ACCOUNT}" --partition=short \
    --array="1-${VALIDATION_TASKS}%${MAX_CONCURRENT}" --job-name=joint_pve_val \
    --output="${VALIDATION_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,JOINT_PVE_SCRIPT_DIR="${VALIDATION_SCRIPT_DIR}",CAL_H2_ENV="${ENV_PATH}",JOINT_PVE_MANIFEST="${VALIDATION_MANIFEST}",JOINT_PVE_CHUNK_MANIFEST="${VALIDATION_CHUNKS}",JOINT_PVE_OUTPUT_ROOT="${VALIDATION_ROOT}",JOINT_PVE_WORK_ROOT="${VALIDATION_ROOT}/work",JOINT_PVE_KEEP_WORK="${KEEP_WORK}" \
    "${VALIDATION_SCRIPT_DIR}/step_joint_pve_chunk.sh")
VALIDATION_JOB_ID=${VALIDATION_JOB%%;*}

FIT_JOB=$(sbatch --parsable --account="${ACCOUNT}" --partition=short \
    --dependency="afterok:${TRAIN_JOB_ID}" --job-name=joint_pve_fit \
    --ntasks=1 --cpus-per-task=1 --mem=24G --time=02:00:00 \
    --output="${TRAIN_ROOT}/logs/%x.%j.log" \
    --wrap="${ENV_PATH}/bin/Rscript ${TRAIN_SCRIPT_DIR}/23_fit_joint_pve_calibrator.R --input_dir=${TRAIN_ROOT}/raw --manifest=${TRAIN_MANIFEST} --config=${TRAIN_ROOT}/config/joint-pve-20260820.tsv --output_dir=${TRAIN_ROOT}/combined")
FIT_JOB_ID=${FIT_JOB%%;*}

EVAL_JOB=$(sbatch --parsable --account="${ACCOUNT}" --partition=short \
    --dependency="afterok:${FIT_JOB_ID}:${VALIDATION_JOB_ID}" \
    --job-name=joint_pve_decide --ntasks=1 --cpus-per-task=1 --mem=24G \
    --time=02:00:00 --output="${VALIDATION_ROOT}/logs/%x.%j.log" \
    --wrap="${ENV_PATH}/bin/Rscript ${VALIDATION_SCRIPT_DIR}/24_evaluate_joint_pve_validation.R --input_dir=${VALIDATION_ROOT}/raw --manifest=${VALIDATION_MANIFEST} --config=${VALIDATION_ROOT}/config/joint-pve-20260820.tsv --criteria=${VALIDATION_ROOT}/config/joint-pve-acceptance-criteria.tsv --model=${TRAIN_ROOT}/combined/joint-pve-calibrator.rds --model_sha256=${TRAIN_ROOT}/combined/joint-pve-calibrator.sha256 --output_dir=${VALIDATION_ROOT}/combined --fail_on_rejection=TRUE")
EVAL_JOB_ID=${EVAL_JOB%%;*}

{
    printf 'stage\tjob_id\n'
    printf 'development_features\t%s\n' "${TRAIN_JOB_ID}"
    printf 'fit_frozen_model\t%s\n' "${FIT_JOB_ID}"
} > "${TRAIN_ROOT}/provenance/submitted-jobs.tsv"
{
    printf 'stage\tjob_id\n'
    printf 'independent_validation_features\t%s\n' "${VALIDATION_JOB_ID}"
    printf 'terminal_evaluation\t%s\n' "${EVAL_JOB_ID}"
} > "${VALIDATION_ROOT}/provenance/submitted-jobs.tsv"

echo "Submitted final joint-PVE experiment"
echo "  development array: ${TRAIN_JOB_ID} (${TRAIN_TASKS} tasks)"
echo "  validation array:  ${VALIDATION_JOB_ID} (${VALIDATION_TASKS} tasks)"
echo "  fit model:          ${FIT_JOB_ID}"
echo "  terminal decision:  ${EVAL_JOB_ID}"
echo "  decision file: ${VALIDATION_ROOT}/combined/joint-pve-terminal-decision.tsv"
