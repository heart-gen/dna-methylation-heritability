#!/bin/bash
# Recover only the terminal evaluator after the original fresh R session did
# not register glmnet's S3 predict method. Simulation outputs and the frozen
# model are read-only inputs; the recovery writes a new immutable derived run.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
TRAIN_ID=${JOINT_PVE_TRAIN_ID:-lgv-joint-pve-train-20260820}
VALIDATION_ID=${JOINT_PVE_VALIDATION_ID:-lgv-joint-pve-validate-20260820}
DECISION_ID=${1:-lgv-joint-pve-decision-20260821}

[[ "${DECISION_ID}" =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "Invalid decision run ID: ${DECISION_ID}" >&2; exit 1;
}
TRAIN_ROOT=${RUN_BASE}/${TRAIN_ID}
VALIDATION_ROOT=${RUN_BASE}/${VALIDATION_ID}
DECISION_ROOT=${RUN_BASE}/${DECISION_ID}
[[ ! -e "${DECISION_ROOT}" ]] || {
    echo "Decision run already exists: ${DECISION_ROOT}" >&2; exit 1;
}

MODEL=${TRAIN_ROOT}/combined/joint-pve-calibrator.rds
MODEL_SHA_FILE=${TRAIN_ROOT}/combined/joint-pve-calibrator.sha256
VALIDATION_RAW=${VALIDATION_ROOT}/raw
VALIDATION_MANIFEST=${VALIDATION_ROOT}/config/scenarios.tsv
CONFIG=${VALIDATION_ROOT}/config/joint-pve-20260820.tsv
CRITERIA=${VALIDATION_ROOT}/config/joint-pve-acceptance-criteria.tsv
for path in "${MODEL}" "${MODEL_SHA_FILE}" "${VALIDATION_MANIFEST}" \
    "${CONFIG}" "${CRITERIA}"; do
    [[ -f "${path}" ]] || { echo "Missing input: ${path}" >&2; exit 1; }
done
EXPECTED=$(($(wc -l < "${VALIDATION_MANIFEST}") - 1))
COMPLETED=$(find "${VALIDATION_RAW}" -maxdepth 1 -type f \
    -name 'scenario-*.tsv' | wc -l)
[[ "${EXPECTED}" -eq 12960 && "${COMPLETED}" -eq "${EXPECTED}" ]] || {
    echo "Validation does not reconcile: expected=${EXPECTED} completed=${COMPLETED}" >&2
    exit 1
}
EXPECTED_SHA=$(tr -d '[:space:]' < "${MODEL_SHA_FILE}")
OBSERVED_SHA=$(sha256sum "${MODEL}" | awk '{print $1}')
[[ "${EXPECTED_SHA}" == "${OBSERVED_SHA}" ]] || {
    echo "Frozen model checksum mismatch" >&2; exit 1;
}

mkdir -p "${DECISION_ROOT}/config" "${DECISION_ROOT}/code" \
    "${DECISION_ROOT}/combined" "${DECISION_ROOT}/logs" \
    "${DECISION_ROOT}/provenance"
cp "${CONFIG}" "${DECISION_ROOT}/config/"
cp "${CRITERIA}" "${DECISION_ROOT}/config/"
cp "${VALIDATION_ROOT}/config/FINAL_JOINT_PVE_STRATEGY.md" \
    "${DECISION_ROOT}/config/"
cp -a "${SCRIPT_DIR}" "${DECISION_ROOT}/code/_h"
RUN_SCRIPT_DIR=${DECISION_ROOT}/code/_h

{
    printf 'field\tvalue\n'
    printf 'run_id\t%s\n' "${DECISION_ID}"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'purpose\ttechnical_recovery_of_terminal_joint_pve_evaluator\n'
    printf 'development_run_id\t%s\n' "${TRAIN_ID}"
    printf 'validation_run_id\t%s\n' "${VALIDATION_ID}"
    printf 'validation_scenarios\t%s\n' "${EXPECTED}"
    printf 'frozen_model_sha256\t%s\n' "${OBSERVED_SHA}"
    printf 'scientific_configuration_changed\tFALSE\n'
    printf 'technical_fix\texplicit_requireNamespace_glmnet_before_S3_predict_dispatch\n'
    printf 'original_failed_job\t9981201\n'
} > "${DECISION_ROOT}/provenance/run-metadata.tsv"
git -C "${REPO_ROOT}" rev-parse HEAD > \
    "${DECISION_ROOT}/provenance/git-commit.txt"

set +e
"${ENV_PATH}/bin/Rscript" \
    "${RUN_SCRIPT_DIR}/24_evaluate_joint_pve_validation.R" \
    --input_dir="${VALIDATION_RAW}" \
    --manifest="${VALIDATION_MANIFEST}" \
    --config="${DECISION_ROOT}/config/joint-pve-20260820.tsv" \
    --criteria="${DECISION_ROOT}/config/joint-pve-acceptance-criteria.tsv" \
    --model="${MODEL}" \
    --model_sha256="${MODEL_SHA_FILE}" \
    --output_dir="${DECISION_ROOT}/combined" \
    --fail_on_rejection=FALSE \
    > "${DECISION_ROOT}/logs/terminal-evaluator.log" 2>&1
STATUS=$?
set -e
printf 'evaluator_exit_status\t%s\n' "${STATUS}" >> \
    "${DECISION_ROOT}/provenance/run-metadata.tsv"
[[ "${STATUS}" -eq 0 ]] || {
    echo "Recovered evaluator failed; inspect ${DECISION_ROOT}/logs/terminal-evaluator.log" >&2
    exit "${STATUS}"
}

find "${DECISION_ROOT}/combined" -maxdepth 1 -type f -print0 | \
    sort -z | xargs -0 sha256sum > \
    "${DECISION_ROOT}/provenance/output-checksums.sha256"
echo "Recovered terminal decision: ${DECISION_ROOT}"
cat "${DECISION_ROOT}/combined/joint-pve-terminal-decision.tsv"
