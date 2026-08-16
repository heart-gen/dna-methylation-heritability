#!/bin/bash

set -euo pipefail

usage() {
    echo "Usage: $0 CALIBRATION_RUN_ID OBSERVED_RUN_ID [POPULATION]"
}
if (( $# < 2 || $# > 3 )); then
    usage >&2
    exit 1
fi

CALIBRATION_RUN_ID=$1
OBSERVED_RUN_ID=$2
POPULATION=${3:-AA}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
MAX_CONCURRENT=${MAX_CONCURRENT:-150}
PLINK_ROOT=${CAL_H2_PLINK_ROOT:-/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/all_individuals}
CALIBRATION_ROOT=${ANALYSIS_DIR}/_m/runs/${CALIBRATION_RUN_ID}
CALIBRATION_MODEL=${CALIBRATION_ROOT}/calibration/elastic-net-calibration.rds
OBSERVED_ROOT=${ANALYSIS_DIR}/_m/observed-runs/${OBSERVED_RUN_ID}
REGIONS=(caudate dlpfc hippocampus)

for id in "${CALIBRATION_RUN_ID}" "${OBSERVED_RUN_ID}"; do
    if [[ ! "${id}" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "Run IDs contain unsupported characters: ${id}" >&2
        exit 1
    fi
done
if [[ "${POPULATION}" != "AA" ]]; then
    echo "Primary calibrated application is restricted to AA; requested ${POPULATION}" >&2
    exit 1
fi
if [[ ! -s "${CALIBRATION_MODEL}" ]]; then
    echo "Calibration model is missing: ${CALIBRATION_MODEL}" >&2
    exit 1
fi
if [[ -e "${OBSERVED_ROOT}" ]]; then
    echo "Observed run already exists: ${OBSERVED_ROOT}" >&2
    exit 1
fi
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

gate_file=$(mktemp /tmp/cal-h2-observed-gate.XXXXXX.tsv)
trap 'rm -f "${gate_file}"' EXIT
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/06_check_acceptance.R" \
    --performance="${CALIBRATION_ROOT}/evaluation/calibration-performance-overall.tsv" \
    --criteria="${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
    --output="${gate_file}" \
    --fail-on-rejection=TRUE

mkdir -p "${OBSERVED_ROOT}/config" "${OBSERVED_ROOT}/code" \
    "${OBSERVED_ROOT}/logs" "${OBSERVED_ROOT}/provenance"
cp -a "${SCRIPT_DIR}" "${OBSERVED_ROOT}/code/_h"
RUN_SCRIPT_DIR=${OBSERVED_ROOT}/code/_h
cp "${CALIBRATION_MODEL}" "${OBSERVED_ROOT}/config/elastic-net-calibration.rds"
cp "${CALIBRATION_ROOT}/evaluation/calibration-performance-overall.tsv" \
    "${OBSERVED_ROOT}/config/calibration-performance-overall.tsv"
cp "${gate_file}" "${OBSERVED_ROOT}/config/calibration-acceptance-results.tsv"
cp "${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
    "${OBSERVED_ROOT}/config/acceptance-criteria.tsv"
CALIBRATION_MODEL=${OBSERVED_ROOT}/config/elastic-net-calibration.rds

printf 'region\tpopulation\texpected_tasks\n' > \
    "${OBSERVED_ROOT}/config/expected-tasks.tsv"
for region in "${REGIONS[@]}"; do
    vmr_file=${REPO_ROOT}/vmr-analysis/all_individuals/${region}/_m/vmr.bed
    if [[ ! -s "${vmr_file}" ]]; then
        echo "VMR file is missing: ${vmr_file}" >&2
        exit 1
    fi
    tasks=$(wc -l < "${vmr_file}")
    printf '%s\t%s\t%s\n' "${region}" "${POPULATION}" "${tasks}" >> \
        "${OBSERVED_ROOT}/config/expected-tasks.tsv"
done
{
    printf 'field\tvalue\n'
    printf 'observed_run_id\t%s\n' "${OBSERVED_RUN_ID}"
    printf 'calibration_run_id\t%s\n' "${CALIBRATION_RUN_ID}"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'population\t%s\n' "${POPULATION}"
    printf 'plink_root\t%s\n' "${PLINK_ROOT}"
    printf 'write_diagnostics\tFALSE\n'
} > "${OBSERVED_ROOT}/provenance/run-metadata.tsv"
git -C "${REPO_ROOT}" rev-parse HEAD > "${OBSERVED_ROOT}/provenance/git-commit.txt"
conda list -p "${ENV_PATH}" --explicit > \
    "${OBSERVED_ROOT}/provenance/conda-explicit-spec.txt"
for path in "${RUN_SCRIPT_DIR}"/*.R "${RUN_SCRIPT_DIR}"/*.sh \
    "${OBSERVED_ROOT}/config"/*; do sha256sum "${path}"; done > \
    "${OBSERVED_ROOT}/provenance/sha256sums.txt"

cd "${REPO_ROOT}"
job_ids=()
printf 'stage\tregion\tjob_id\tdependency\n' > \
    "${OBSERVED_ROOT}/provenance/submitted-jobs.tsv"
while IFS=$'\t' read -r region population tasks; do
    [[ "${region}" == "region" ]] && continue
    job=$(sbatch --parsable --account="${ACCOUNT}" \
        --array="1-${tasks}%${MAX_CONCURRENT}" \
        --job-name="cal_h2_${region}" \
        --output="${OBSERVED_ROOT}/logs/%x.%A_%a.log" \
        --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_REPO_ROOT="${REPO_ROOT}",CAL_H2_CALIBRATION_MODEL="${CALIBRATION_MODEL}",CAL_H2_OBSERVED_OUTPUT_ROOT="${OBSERVED_ROOT}/results",CAL_H2_PLINK_ROOT="${PLINK_ROOT}",CAL_H2_WRITE_DIAGNOSTICS=FALSE,REGION="${region}",POPULATION="${population}" \
        "${RUN_SCRIPT_DIR}/step_5_estimate_observed_vmr.sh")
    job_id=${job%%;*}
    job_ids+=("${job_id}")
    printf 'estimate\t%s\t%s\tNA\n' "${region}" "${job_id}" >> \
        "${OBSERVED_ROOT}/provenance/submitted-jobs.tsv"
done < "${OBSERVED_ROOT}/config/expected-tasks.tsv"
dependency=$(IFS=:; echo "${job_ids[*]}")
qc_job=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${dependency}" \
    --output="${OBSERVED_ROOT}/logs/%x.%j.log" \
    --export=ALL,CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_OBSERVED_OUTPUT_ROOT="${OBSERVED_ROOT}/results",CAL_H2_EXPECTED_TASKS="${OBSERVED_ROOT}/config/expected-tasks.tsv" \
    "${RUN_SCRIPT_DIR}/step_6_combine_observed.sh")
qc_job_id=${qc_job%%;*}
printf 'combine_qc\tall\t%s\tafterok:%s\n' "${qc_job_id}" "${dependency}" >> \
    "${OBSERVED_ROOT}/provenance/submitted-jobs.tsv"

echo "Observed run: ${OBSERVED_ROOT}"
echo "Region arrays: ${job_ids[*]}"
echo "Combine/QC job: ${qc_job_id}"
