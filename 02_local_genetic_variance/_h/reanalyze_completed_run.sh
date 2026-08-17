#!/bin/bash

set -euo pipefail

usage() {
    echo "Usage: $0 SOURCE_RUN_ID DERIVED_RUN_ID"
}
if (( $# != 2 )); then
    usage >&2
    exit 1
fi

SOURCE_RUN_ID=$1
DERIVED_RUN_ID=$2
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
RUN_BASE=${CAL_H2_RUN_BASE:-${ANALYSIS_DIR}/_m/runs}
SOURCE_ROOT=${RUN_BASE}/${SOURCE_RUN_ID}
DERIVED_ROOT=${RUN_BASE}/${DERIVED_RUN_ID}

for id in "${SOURCE_RUN_ID}" "${DERIVED_RUN_ID}"; do
    if [[ ! "${id}" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "Run IDs may contain only letters, numbers, periods, underscores, and hyphens" >&2
        exit 1
    fi
done
if [[ ! -s "${SOURCE_ROOT}/config/scenarios.tsv" ]]; then
    echo "Source scenario manifest is missing: ${SOURCE_ROOT}/config/scenarios.tsv" >&2
    exit 1
fi
if [[ -e "${DERIVED_ROOT}" ]]; then
    echo "Derived run already exists: ${DERIVED_ROOT}" >&2
    exit 1
fi
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

expected_calibration=$(awk -F '\t' 'NR > 1 && $2 == "calibration" {n++} END {print n+0}' \
    "${SOURCE_ROOT}/config/scenarios.tsv")
expected_evaluation=$(awk -F '\t' 'NR > 1 && $2 == "evaluation" {n++} END {print n+0}' \
    "${SOURCE_ROOT}/config/scenarios.tsv")
actual_calibration=$(find "${SOURCE_ROOT}/raw/calibration" -maxdepth 1 \
    -type f -name 'scenario-*.tsv' | wc -l)
actual_evaluation=$(find "${SOURCE_ROOT}/raw/evaluation" -maxdepth 1 \
    -type f -name 'scenario-*.tsv' | wc -l)
if (( actual_calibration != expected_calibration ||
      actual_evaluation != expected_evaluation )); then
    echo "Source raw results are incomplete" >&2
    echo "Calibration: ${actual_calibration}/${expected_calibration}" >&2
    echo "Evaluation: ${actual_evaluation}/${expected_evaluation}" >&2
    exit 1
fi

mkdir -p "${DERIVED_ROOT}/config" "${DERIVED_ROOT}/code" \
    "${DERIVED_ROOT}/provenance" "${DERIVED_ROOT}/logs"
cp -a "${SCRIPT_DIR}" "${DERIVED_ROOT}/code/_h"
RUN_SCRIPT_DIR=${DERIVED_ROOT}/code/_h
cp "${SOURCE_ROOT}/config/analysis.tsv" \
    "${DERIVED_ROOT}/config/source-analysis.tsv"
cp "${SOURCE_ROOT}/config/scenarios.tsv" \
    "${DERIVED_ROOT}/config/source-scenarios.tsv"
cp "${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
    "${DERIVED_ROOT}/config/acceptance-criteria.tsv"
cp "${ANALYSIS_DIR}/environment.yml" "${DERIVED_ROOT}/config/environment.yml"
ln -s "${SOURCE_ROOT}/raw" "${DERIVED_ROOT}/raw"

{
    printf 'field\tvalue\n'
    printf 'source_run\t%s\n' "${SOURCE_RUN_ID}"
    printf 'derived_run\t%s\n' "${DERIVED_RUN_ID}"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'reason\tnull_aware_forward_hybrid_calibration\n'
    printf 'null_alpha\t0.05\n'
    printf 'null_threshold_method\tsplit_conformal_order_statistic\n'
    printf 'source_calibration_results\t%s\n' "${actual_calibration}"
    printf 'source_evaluation_results\t%s\n' "${actual_evaluation}"
} > "${DERIVED_ROOT}/provenance/derivation.tsv"
git -C "${REPO_ROOT}" rev-parse HEAD > "${DERIVED_ROOT}/provenance/git-commit.txt"
conda list -p "${ENV_PATH}" --explicit > \
    "${DERIVED_ROOT}/provenance/conda-explicit-spec.txt"
for path in "${RUN_SCRIPT_DIR}"/*.R "${RUN_SCRIPT_DIR}"/*.sh \
    "${DERIVED_ROOT}/config"/*; do
    sha256sum "${path}"
done > "${DERIVED_ROOT}/provenance/sha256sums.txt"

"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/02_fit_calibration.R" \
    --input="${DERIVED_ROOT}/raw/calibration" \
    --criteria="${DERIVED_ROOT}/config/acceptance-criteria.tsv" \
    --output-model="${DERIVED_ROOT}/calibration/elastic-net-calibration.rds" \
    --output-manifest="${DERIVED_ROOT}/calibration/calibration-manifest.tsv" \
    --output-tuning="${DERIVED_ROOT}/calibration/hybrid-weight-tuning.tsv" \
    --session-info="${DERIVED_ROOT}/calibration/session-info.txt" \
    > "${DERIVED_ROOT}/logs/fit-calibration.log" 2>&1
"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/03_evaluate_calibration.R" \
    --input="${DERIVED_ROOT}/raw/evaluation" \
    --model="${DERIVED_ROOT}/calibration/elastic-net-calibration.rds" \
    --output-dir="${DERIVED_ROOT}/evaluation" \
    > "${DERIVED_ROOT}/logs/evaluate-calibration.log" 2>&1
"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/05_plot_calibration.R" \
    --input="${DERIVED_ROOT}/evaluation/calibrated-evaluation-estimates.tsv" \
    --output-dir="${DERIVED_ROOT}/figures" \
    > "${DERIVED_ROOT}/logs/plot-calibration.log" 2>&1
"${ENV_PATH}/bin/Rscript" "${RUN_SCRIPT_DIR}/06_check_acceptance.R" \
    --performance="${DERIVED_ROOT}/evaluation/calibration-performance-overall.tsv" \
    --criteria="${DERIVED_ROOT}/config/acceptance-criteria.tsv" \
    --model="${DERIVED_ROOT}/calibration/elastic-net-calibration.rds" \
    --output="${DERIVED_ROOT}/evaluation/acceptance-results.tsv" \
    --fail-on-rejection=TRUE \
    > "${DERIVED_ROOT}/logs/check-acceptance.log" 2>&1

echo "Derived calibration run completed: ${DERIVED_ROOT}"
