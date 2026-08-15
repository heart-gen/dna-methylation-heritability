#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
RSCRIPT=${ENV_PATH}/bin/Rscript
SMOKE_DIR=$(mktemp -d /tmp/calibrated-local-h2-smoke.XXXXXX)
trap 'rm -rf "${SMOKE_DIR}"' EXIT

"${RSCRIPT}" "${SCRIPT_DIR}/test_functions.R"
"${RSCRIPT}" "${ANALYSIS_DIR}/_h/00_make_manifest.R" \
    --config="${ANALYSIS_DIR}/config/smoke.tsv" \
    --output="${SMOKE_DIR}/scenarios.tsv"
TASKS=$(($(wc -l < "${SMOKE_DIR}/scenarios.tsv") - 1))
for TASK_ID in $(seq 1 "${TASKS}"); do
    "${RSCRIPT}" "${ANALYSIS_DIR}/_h/01_simulate_and_crossfit.R" \
        --manifest="${SMOKE_DIR}/scenarios.tsv" \
        --task-id="${TASK_ID}" \
        --output-root="${SMOKE_DIR}/raw"
done
"${RSCRIPT}" "${ANALYSIS_DIR}/_h/02_fit_calibration.R" \
    --input="${SMOKE_DIR}/raw/calibration" \
    --output-model="${SMOKE_DIR}/calibration.rds" \
    --output-manifest="${SMOKE_DIR}/calibration-manifest.tsv" \
    --session-info="${SMOKE_DIR}/calibration-session.txt"
"${RSCRIPT}" "${ANALYSIS_DIR}/_h/03_evaluate_calibration.R" \
    --input="${SMOKE_DIR}/raw/evaluation" \
    --model="${SMOKE_DIR}/calibration.rds" \
    --output-dir="${SMOKE_DIR}/evaluation"
"${RSCRIPT}" "${ANALYSIS_DIR}/_h/05_plot_calibration.R" \
    --input="${SMOKE_DIR}/evaluation/calibrated-evaluation-estimates.tsv" \
    --output-dir="${SMOKE_DIR}/figures"
"${RSCRIPT}" "${ANALYSIS_DIR}/_h/06_check_acceptance.R" \
    --performance="${SMOKE_DIR}/evaluation/calibration-performance-overall.tsv" \
    --criteria="${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
    --output="${SMOKE_DIR}/evaluation/acceptance-results.tsv"
test -s "${SMOKE_DIR}/evaluation/calibrated-evaluation-estimates.tsv"
test -s "${SMOKE_DIR}/evaluation/calibration-performance-overall.tsv"
test -s "${SMOKE_DIR}/evaluation/acceptance-results.tsv"
test -s "${SMOKE_DIR}/figures/calibration-truth-versus-estimate.pdf"

# Exercise the final observed-data combiner with one complete, in-domain VMR
# result from each prespecified brain region.
OBSERVED_DIR=${SMOKE_DIR}/observed
EXPECTED_FILE=${SMOKE_DIR}/expected-observed-tasks.tsv
printf 'region\tpopulation\texpected_tasks\n' > "${EXPECTED_FILE}"
for REGION in caudate dlpfc hippocampus; do
    mkdir -p "${OBSERVED_DIR}/${REGION}/AA/summary"
    printf '%s\tAA\t1\n' "${REGION}" >> "${EXPECTED_FILE}"
    {
        printf 'task_id\tregion\tpopulation\tchromosome\tstart\tcalibration_status\tpositive_signal\n'
        printf '1\t%s\tAA\tchr1\t100\twithin_domain\tTRUE\n' "${REGION}"
    } > "${OBSERVED_DIR}/${REGION}/AA/summary/vmr-0000001.tsv"
done
"${RSCRIPT}" "${ANALYSIS_DIR}/_h/07_combine_observed.R" \
    --input="${OBSERVED_DIR}" \
    --expected="${EXPECTED_FILE}" \
    --output-dir="${SMOKE_DIR}/observed-combined"
test "$(($(wc -l < "${SMOKE_DIR}/observed-combined/calibrated-local-h2-AA-vmrs.tsv") - 1))" -eq 3
test "$(($(wc -l < "${SMOKE_DIR}/observed-combined/observed-run-qc.tsv") - 1))" -eq 3
"${RSCRIPT}" -e \
    'x <- read.delim(commandArgs(TRUE)[1]); stopifnot(all(x$overall_qc_pass))' \
    "${SMOKE_DIR}/observed-combined/observed-run-qc.tsv"

# A computational failure is never absorbed by the locus-QC tolerance.
FAILED_OBSERVED_DIR=${SMOKE_DIR}/observed-with-computational-failure
FAILED_EXPECTED_FILE=${SMOKE_DIR}/expected-failed-observed-tasks.tsv
mkdir -p "${FAILED_OBSERVED_DIR}/caudate/AA/summary" \
    "${FAILED_OBSERVED_DIR}/caudate/AA/failures"
printf 'region\tpopulation\texpected_tasks\ncaudate\tAA\t2\n' > \
    "${FAILED_EXPECTED_FILE}"
{
    printf 'task_id\tregion\tpopulation\tchromosome\tstart\tcalibration_status\tpositive_signal\n'
    printf '1\tcaudate\tAA\tchr1\t100\twithin_domain\tTRUE\n'
} > "${FAILED_OBSERVED_DIR}/caudate/AA/summary/vmr-0000001.tsv"
{
    printf 'task_id\tregion\tpopulation\texit_status\tlog_file\n'
    printf '2\tcaudate\tAA\t1\ttest.log\n'
} > "${FAILED_OBSERVED_DIR}/caudate/AA/failures/vmr-0000002.tsv"
"${RSCRIPT}" "${ANALYSIS_DIR}/_h/07_combine_observed.R" \
    --input="${FAILED_OBSERVED_DIR}" \
    --expected="${FAILED_EXPECTED_FILE}" \
    --output-dir="${SMOKE_DIR}/observed-failure-combined" \
    --fail-on-qc=FALSE
"${RSCRIPT}" -e \
    'x <- read.delim(commandArgs(TRUE)[1]); stopifnot(x$failure_qc_pass, !x$computational_qc_pass, !x$overall_qc_pass, x$qc_failed_tasks == 0, x$computational_failed_tasks == 1)' \
    "${SMOKE_DIR}/observed-failure-combined/observed-run-qc.tsv"
echo "End-to-end smoke test passed"
