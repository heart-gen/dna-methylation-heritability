#!/bin/bash
# Resume scheduler submission after the immutable final joint-PVE run roots
# were prepared but no scheduler job was created.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
TRAIN_ID=${1:-lgv-joint-pve-train-20260820}
VALIDATION_ID=${2:-lgv-joint-pve-validate-20260820}
TRAIN_ROOT=${RUN_BASE}/${TRAIN_ID}
VALIDATION_ROOT=${RUN_BASE}/${VALIDATION_ID}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
KEEP_WORK=${JOINT_PVE_KEEP_WORK:-FALSE}

for root in "${TRAIN_ROOT}" "${VALIDATION_ROOT}"; do
    [[ -d "${root}" ]] || { echo "Missing prepared run: ${root}" >&2; exit 1; }
    [[ ! -e "${root}/provenance/submitted-jobs.tsv" ]] || {
        echo "Submission metadata already exists: ${root}" >&2; exit 1;
    }
done
if find "${TRAIN_ROOT}/raw" "${VALIDATION_ROOT}/raw" -type f -print -quit |
    grep -q .; then
    echo "Prepared runs already contain raw outputs; refusing ambiguous resume" >&2
    exit 1
fi

TRAIN_SCRIPT_DIR=${TRAIN_ROOT}/code/_h
VALIDATION_SCRIPT_DIR=${VALIDATION_ROOT}/code/_h
TRAIN_MANIFEST=${TRAIN_ROOT}/config/scenarios.tsv
VALIDATION_MANIFEST=${VALIDATION_ROOT}/config/scenarios.tsv
TRAIN_CHUNKS=${TRAIN_ROOT}/config/chunk-manifest.tsv
VALIDATION_CHUNKS=${VALIDATION_ROOT}/config/chunk-manifest.tsv
CONFIG=${TRAIN_ROOT}/config/joint-pve-20260820.tsv
for file in "${TRAIN_MANIFEST}" "${VALIDATION_MANIFEST}" "${TRAIN_CHUNKS}" \
    "${VALIDATION_CHUNKS}" "${CONFIG}"; do
    [[ -f "${file}" ]] || { echo "Missing prepared file: ${file}" >&2; exit 1; }
done
config_value() {
    awk -F '\t' -v key="$1" '$1==key {print $2; exit}' "${CONFIG}"
}
DEVELOPMENT_REL=$(config_value development_bslmm_run_relpath)
DEVELOPMENT_RAW=${ANALYSIS_DIR}/${DEVELOPMENT_REL}/raw
[[ -d "${DEVELOPMENT_RAW}" ]] || {
    echo "Missing development BSLMM raw directory" >&2; exit 1;
}
MAX_CONCURRENT=${MAX_CONCURRENT:-$(config_value max_concurrent)}
TRAIN_TASKS=$(awk -F '\t' 'NR>1 {last=$1} END {print last+0}' "${TRAIN_CHUNKS}")
VALIDATION_TASKS=$(awk -F '\t' 'NR>1 {last=$1} END {print last+0}' "${VALIDATION_CHUNKS}")
[[ "${TRAIN_TASKS}" -gt 0 && "${VALIDATION_TASKS}" -gt 0 ]] || {
    echo "Prepared chunk manifests are empty" >&2; exit 1;
}

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

echo "Resumed final joint-PVE submission"
echo "  development array: ${TRAIN_JOB_ID}"
echo "  validation array:  ${VALIDATION_JOB_ID}"
echo "  fit model:          ${FIT_JOB_ID}"
echo "  terminal decision:  ${EVAL_JOB_ID}"
