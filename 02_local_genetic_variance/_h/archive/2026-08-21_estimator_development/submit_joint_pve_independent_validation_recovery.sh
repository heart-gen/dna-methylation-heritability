#!/bin/bash
# Submit the one seed-independence recovery validation against the already
# frozen final joint-PVE model. No model fitting or gate modification occurs.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
KEEP_WORK=${JOINT_PVE_KEEP_WORK:-FALSE}
DRY_RUN=${SUBMIT_JOINT_PVE_DRY_RUN:-FALSE}

BASE_CONFIG=${ANALYSIS_DIR}/config/joint-pve-20260820.tsv
CRITERIA=${ANALYSIS_DIR}/config/joint-pve-acceptance-criteria.tsv
STRATEGY=${ANALYSIS_DIR}/config/FINAL_JOINT_PVE_STRATEGY.md
RECOVERY_LOCK=${ANALYSIS_DIR}/config/joint-pve-validation-recovery-20260821.tsv
RECOVERY_NOTE=${ANALYSIS_DIR}/config/JOINT_PVE_VALIDATION_INDEPENDENCE_RECOVERY.md
TRAIN_ID=lgv-joint-pve-train-20260820
TRAIN_ROOT=${JOINT_PVE_TRAIN_ROOT:-${ANALYSIS_DIR}/_m/runs/${TRAIN_ID}}
MODEL=${TRAIN_ROOT}/combined/joint-pve-calibrator.rds
MODEL_SHA_FILE=${TRAIN_ROOT}/combined/joint-pve-calibrator.sha256

lock_value() {
    awk -F '\t' -v key="$1" '$1==key {print $2; exit}' "${RECOVERY_LOCK}"
}
VALIDATION_ID=$(lock_value validation_run_id)
DECISION_ID=$(lock_value decision_run_id)
VALIDATION_ROOT=${RUN_BASE}/${VALIDATION_ID}
DECISION_ROOT=${RUN_BASE}/${DECISION_ID}

for path in "${BASE_CONFIG}" "${CRITERIA}" "${STRATEGY}" \
    "${RECOVERY_LOCK}" "${RECOVERY_NOTE}" "${MODEL}" "${MODEL_SHA_FILE}"; do
    [[ -f "${path}" ]] || { echo "Missing locked input: ${path}" >&2; exit 1; }
done
export PATH="/projects/p32505/opt/bin:${PATH}"
command -v gemma >/dev/null || { echo "gemma not on PATH" >&2; exit 1; }
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"
for root in "${VALIDATION_ROOT}" "${DECISION_ROOT}"; do
    [[ ! -e "${root}" ]] || { echo "Run already exists: ${root}" >&2; exit 1; }
done

check_sha() {
    local path=$1
    local expected=$2
    local observed
    observed=$(sha256sum "${path}" | awk '{print $1}')
    [[ "${observed}" == "${expected}" ]] || {
        echo "Checksum mismatch for ${path}" >&2; exit 1;
    }
}
check_sha "${BASE_CONFIG}" "$(lock_value base_config_sha256)"
check_sha "${CRITERIA}" "$(lock_value criteria_sha256)"
check_sha "${MODEL}" "$(lock_value frozen_model_sha256)"
[[ "$(tr -d '[:space:]' < "${MODEL_SHA_FILE}")" == \
   "$(lock_value frozen_model_sha256)" ]] || {
    echo "Frozen model checksum file mismatch" >&2; exit 1;
}

TMP_DIR=$(mktemp -d /tmp/joint-pve-validation-recovery.XXXXXX)
trap 'rm -rf "${TMP_DIR}"' EXIT
EFFECTIVE_CONFIG=${TMP_DIR}/joint-pve-effective-recovery.tsv
CANDIDATE_MANIFEST=${TMP_DIR}/scenarios.tsv
PREFLIGHT_SEED_AUDIT=${TMP_DIR}/seed-independence-audit.tsv

"${ENV_PATH}/bin/Rscript" \
    "${SCRIPT_DIR}/25_prepare_joint_pve_validation_recovery.R" \
    --base_config="${BASE_CONFIG}" --recovery_lock="${RECOVERY_LOCK}" \
    --output="${EFFECTIVE_CONFIG}"
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/21_make_joint_pve_manifest.R" \
    --config="${EFFECTIVE_CONFIG}" --mode=validation \
    --output="${CANDIDATE_MANIFEST}"
"${ENV_PATH}/bin/Rscript" \
    "${SCRIPT_DIR}/26_audit_joint_pve_seed_independence.R" \
    --candidate="${CANDIDATE_MANIFEST}" --repo_root="${REPO_ROOT}" \
    --output="${PREFLIGHT_SEED_AUDIT}"

EXPECTED=$(($(wc -l < "${CANDIDATE_MANIFEST}") - 1))
[[ "${EXPECTED}" -eq 12960 ]] || {
    echo "Expected 12960 validation scenarios, found ${EXPECTED}" >&2; exit 1;
}
SCENARIOS_PER_TASK=$(awk -F '\t' '$1=="scenarios_per_array_task" {print $2; exit}' "${EFFECTIVE_CONFIG}")
MAX_CONCURRENT=$(awk -F '\t' '$1=="max_concurrent" {print $2; exit}' "${EFFECTIVE_CONFIG}")
TASKS=$(( (EXPECTED + SCENARIOS_PER_TASK - 1) / SCENARIOS_PER_TASK ))

for root in "${VALIDATION_ROOT}" "${DECISION_ROOT}"; do
    mkdir -p "${root}/config" "${root}/logs" "${root}/provenance" \
        "${root}/code" "${root}/combined"
    cp "${BASE_CONFIG}" "${root}/config/"
    cp "${EFFECTIVE_CONFIG}" "${root}/config/"
    cp "${CRITERIA}" "${root}/config/"
    cp "${STRATEGY}" "${root}/config/"
    cp "${RECOVERY_LOCK}" "${root}/config/"
    cp "${RECOVERY_NOTE}" "${root}/config/"
    cp -a "${SCRIPT_DIR}" "${root}/code/_h"
done
mkdir -p "${VALIDATION_ROOT}/raw" "${VALIDATION_ROOT}/work"
cp "${CANDIDATE_MANIFEST}" "${VALIDATION_ROOT}/config/scenarios.tsv"
"${ENV_PATH}/bin/Rscript" \
    "${VALIDATION_ROOT}/code/_h/26_audit_joint_pve_seed_independence.R" \
    --candidate="${VALIDATION_ROOT}/config/scenarios.tsv" \
    --repo_root="${REPO_ROOT}" \
    --output="${VALIDATION_ROOT}/config/seed-independence-audit.tsv"
cp "${VALIDATION_ROOT}/config/seed-independence-audit.tsv" \
    "${DECISION_ROOT}/config/"

CHUNK_MANIFEST=${VALIDATION_ROOT}/config/chunk-manifest.tsv
{
    printf 'chunk_id\tscenario_id\n'
    awk -F '\t' -v size="${SCENARIOS_PER_TASK}" \
        'NR > 1 {print int((NR - 2) / size) + 1 "\t" $1}' \
        "${VALIDATION_ROOT}/config/scenarios.tsv"
} > "${CHUNK_MANIFEST}"

COMMIT=$(git -C "${REPO_ROOT}" rev-parse HEAD)
for root in "${VALIDATION_ROOT}" "${DECISION_ROOT}"; do
    printf '%s\n' "${COMMIT}" > "${root}/provenance/git-commit.txt"
done
{
    printf 'field\tvalue\n'
    printf 'run_id\t%s\n' "${VALIDATION_ID}"
    printf 'purpose\tfully_independent_frozen_joint_pve_validation\n'
    printf 'n_scenarios\t%s\n' "${EXPECTED}"
    printf 'prior_validation_status\tinvalid_for_acceptance_seed_overlap\n'
    printf 'prior_overlap_seed_count\t384\n'
    printf 'frozen_training_run_id\t%s\n' "${TRAIN_ID}"
    printf 'frozen_model_sha256\t%s\n' "$(lock_value frozen_model_sha256)"
    printf 'criteria_sha256\t%s\n' "$(lock_value criteria_sha256)"
    printf 'effective_config_sha256\t%s\n' "$(sha256sum "${EFFECTIVE_CONFIG}" | awk '{print $1}')"
    printf 'manifest_sha256\t%s\n' "$(sha256sum "${CANDIDATE_MANIFEST}" | awk '{print $1}')"
    printf 'seed_overlap_with_prior_manifests\t0\n'
    printf 'model_changed\tFALSE\n'
    printf 'criteria_changed\tFALSE\n'
} > "${VALIDATION_ROOT}/provenance/run-metadata.tsv"
{
    printf 'field\tvalue\n'
    printf 'run_id\t%s\n' "${DECISION_ID}"
    printf 'purpose\tterminal_decision_from_independent_recovery_validation\n'
    printf 'validation_run_id\t%s\n' "${VALIDATION_ID}"
    printf 'frozen_training_run_id\t%s\n' "${TRAIN_ID}"
    printf 'frozen_model_sha256\t%s\n' "$(lock_value frozen_model_sha256)"
    printf 'criteria_sha256\t%s\n' "$(lock_value criteria_sha256)"
    printf 'model_changed\tFALSE\n'
    printf 'criteria_changed\tFALSE\n'
} > "${DECISION_ROOT}/provenance/run-metadata.tsv"

if [[ "${DRY_RUN,,}" == "true" ]]; then
    echo "Prepared independent validation recovery (dry run)"
    echo "  validation: ${VALIDATION_ID} (${EXPECTED} scenarios / ${TASKS} tasks)"
    echo "  decision:   ${DECISION_ID}"
    exit 0
fi

VALIDATION_SCRIPT_DIR=${VALIDATION_ROOT}/code/_h
DECISION_SCRIPT_DIR=${DECISION_ROOT}/code/_h
VALIDATION_JOB=$(sbatch --parsable --account="${ACCOUNT}" --partition=short \
    --array="1-${TASKS}%${MAX_CONCURRENT}" --job-name=joint_pve_val2 \
    --output="${VALIDATION_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,JOINT_PVE_SCRIPT_DIR="${VALIDATION_SCRIPT_DIR}",CAL_H2_ENV="${ENV_PATH}",JOINT_PVE_MANIFEST="${VALIDATION_ROOT}/config/scenarios.tsv",JOINT_PVE_CHUNK_MANIFEST="${CHUNK_MANIFEST}",JOINT_PVE_OUTPUT_ROOT="${VALIDATION_ROOT}",JOINT_PVE_WORK_ROOT="${VALIDATION_ROOT}/work",JOINT_PVE_KEEP_WORK="${KEEP_WORK}" \
    "${VALIDATION_SCRIPT_DIR}/step_joint_pve_chunk.sh")
VALIDATION_JOB_ID=${VALIDATION_JOB%%;*}

DECISION_JOB=$(sbatch --parsable --account="${ACCOUNT}" --partition=short \
    --dependency="afterok:${VALIDATION_JOB_ID}" --job-name=joint_pve_decide2 \
    --output="${DECISION_ROOT}/logs/%x.%j.log" \
    --export=ALL,JOINT_PVE_SCRIPT_DIR="${DECISION_SCRIPT_DIR}",CAL_H2_ENV="${ENV_PATH}",JOINT_PVE_VALIDATION_ROOT="${VALIDATION_ROOT}",JOINT_PVE_DECISION_ROOT="${DECISION_ROOT}",JOINT_PVE_TRAIN_ROOT="${TRAIN_ROOT}" \
    "${DECISION_SCRIPT_DIR}/step_joint_pve_terminal_evaluation.sh")
DECISION_JOB_ID=${DECISION_JOB%%;*}

{
    printf 'stage\tjob_id\n'
    printf 'independent_validation_features\t%s\n' "${VALIDATION_JOB_ID}"
} > "${VALIDATION_ROOT}/provenance/submitted-jobs.tsv"
{
    printf 'stage\tjob_id\n'
    printf 'terminal_evaluation\t%s\n' "${DECISION_JOB_ID}"
} > "${DECISION_ROOT}/provenance/submitted-jobs.tsv"

echo "Submitted independent validation recovery"
echo "  validation array: ${VALIDATION_JOB_ID} (${TASKS} tasks)"
echo "  terminal decision: ${DECISION_JOB_ID}"
echo "  decision run: ${DECISION_ROOT}"
