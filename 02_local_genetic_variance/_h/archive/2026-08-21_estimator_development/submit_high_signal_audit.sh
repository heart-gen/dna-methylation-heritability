#!/bin/bash

set -euo pipefail
if (( $# != 2 )); then
    echo "Usage: $0 RUN_ID INITIAL_COMBINED_OBSERVED_TSV" >&2
    exit 1
fi
RUN_ID=$1
INPUT=$2
[[ "${RUN_ID}" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Invalid run ID" >&2; exit 1; }
[[ -s "${INPUT}" ]] || { echo "Missing observed input: ${INPUT}" >&2; exit 1; }

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
MAX_CONCURRENT=${MAX_CONCURRENT:-100}
MODEL=${CAL_H2_CALIBRATION_MODEL:-${ANALYSIS_DIR}/_m/calibration_frozen/elastic-net-calibration.rds}
RUN_ROOT=${ANALYSIS_DIR}/_m/runs/${RUN_ID}
[[ ! -e "${RUN_ROOT}" ]] || { echo "Run already exists: ${RUN_ROOT}" >&2; exit 1; }
[[ -s "${MODEL}" ]] || { echo "Missing calibration model: ${MODEL}" >&2; exit 1; }

mkdir -p "${RUN_ROOT}"/{code,config,logs,provenance,summaries,tasks}
cp -a "${SCRIPT_DIR}" "${RUN_ROOT}/code/_h"
RUN_SCRIPT_DIR=${RUN_ROOT}/code/_h
cp "${INPUT}" "${RUN_ROOT}/config/initial-observed-combined.tsv"
cp "${MODEL}" "${RUN_ROOT}/config/elastic-net-calibration.rds"
"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/08_make_high_signal_audit_manifest.R" \
    --input="${RUN_ROOT}/config/initial-observed-combined.tsv" \
    --output="${RUN_ROOT}/config/high-signal-audit-manifest.tsv"
TASKS=$(($(wc -l < "${RUN_ROOT}/config/high-signal-audit-manifest.tsv") - 1))

git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_ROOT}/provenance/git-commit.txt"
sha256sum "${RUN_ROOT}/config/"* > "${RUN_ROOT}/provenance/input-sha256sums.txt"
if [[ -n "${KINSHIP_TABLE:-}" && -s "${KINSHIP_TABLE}" ]]; then
    cp "${KINSHIP_TABLE}" "${RUN_ROOT}/provenance/king-related-pairs.kin0"
fi

AUDIT_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --array="1-${TASKS}%${MAX_CONCURRENT}" \
    --output="${RUN_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_REPO_ROOT="${REPO_ROOT}",CAL_H2_AUDIT_MANIFEST="${RUN_ROOT}/config/high-signal-audit-manifest.tsv",CAL_H2_AUDIT_OUTPUT_ROOT="${RUN_ROOT}",CAL_H2_CALIBRATION_MODEL="${RUN_ROOT}/config/elastic-net-calibration.rds" \
    "${RUN_SCRIPT_DIR}/step_7_high_signal_audit.sh")
AUDIT_JOB_ID=${AUDIT_JOB%%;*}
SUMMARY_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${AUDIT_JOB_ID}" \
    --output="${RUN_ROOT}/logs/%x.%j.log" \
    --export=ALL,CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_AUDIT_RUN_ROOT="${RUN_ROOT}" \
    "${RUN_SCRIPT_DIR}/step_8_summarize_high_signal_audit.sh")
SUMMARY_JOB_ID=${SUMMARY_JOB%%;*}
printf 'stage\tjob_id\tdependency\nfold_relatedness_audit\t%s\tNA\nsummarize\t%s\tafterok:%s\n' \
    "${AUDIT_JOB_ID}" "${SUMMARY_JOB_ID}" "${AUDIT_JOB_ID}" > \
    "${RUN_ROOT}/provenance/submitted-jobs.tsv"
echo "Audit run: ${RUN_ROOT}"
echo "Audit array: ${AUDIT_JOB_ID}; summary: ${SUMMARY_JOB_ID}"
