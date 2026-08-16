#!/bin/bash

set -euo pipefail

usage() {
    echo "Usage: $0 SOURCE_OBSERVED_RUN_ID RECOVERY_RUN_ID"
    echo "Reruns only source computational-failure records under corrected classification."
}
if (( $# != 2 )); then
    usage >&2
    exit 1
fi

SOURCE_RUN_ID=$1
RECOVERY_RUN_ID=$2
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
MAX_CONCURRENT=${MAX_CONCURRENT:-150}
RUN_BASE=${CAL_H2_OBSERVED_RUN_BASE:-${ANALYSIS_DIR}/_m/observed-runs}
SOURCE_ROOT=${RUN_BASE}/${SOURCE_RUN_ID}
RECOVERY_ROOT=${RUN_BASE}/${RECOVERY_RUN_ID}
REGIONS=(caudate dlpfc hippocampus)

for id in "${SOURCE_RUN_ID}" "${RECOVERY_RUN_ID}"; do
    if [[ ! "${id}" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "Run IDs contain unsupported characters: ${id}" >&2
        exit 1
    fi
done
if [[ ! -s "${SOURCE_ROOT}/config/expected-tasks.tsv" ||
      ! -s "${SOURCE_ROOT}/config/elastic-net-calibration.rds" ]]; then
    echo "Source observed run is missing required configuration" >&2
    exit 1
fi
if [[ -e "${RECOVERY_ROOT}" ]]; then
    echo "Recovery run already exists: ${RECOVERY_ROOT}" >&2
    exit 1
fi
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

SOURCE_METADATA=${SOURCE_ROOT}/provenance/run-metadata.tsv
if [[ ! -s "${SOURCE_METADATA}" ]]; then
    SOURCE_METADATA=${SOURCE_ROOT}/provenance/recovery-metadata.tsv
fi
if [[ ! -s "${SOURCE_METADATA}" ]]; then
    echo "Source observed run has no run or recovery metadata" >&2
    exit 1
fi
POPULATION=$(awk -F '\t' '$1 == "population" {print $2}' "${SOURCE_METADATA}")
PLINK_ROOT=$(awk -F '\t' '$1 == "plink_root" {print $2}' "${SOURCE_METADATA}")
if [[ "${POPULATION}" != "AA" || -z "${PLINK_ROOT}" ]]; then
    echo "Recovery currently requires an AA source run with a recorded PLINK root" >&2
    exit 1
fi
PHENOTYPE_ROOT=${CAL_H2_PHENOTYPE_ROOT:-/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/all_individuals}
PGEN_ROOT=${CAL_H2_PGEN_ROOT:-/projects/b1213/users/alexis/projects/dna-methylation-heritability/inputs/genotypes/all_individuals}
RECOVERED_PLINK_ROOT=${RECOVERY_ROOT}/recovered-inputs

mkdir -p "${RECOVERY_ROOT}/config" "${RECOVERY_ROOT}/code" \
    "${RECOVERY_ROOT}/logs" "${RECOVERY_ROOT}/provenance"
cp -a "${SCRIPT_DIR}" "${RECOVERY_ROOT}/code/_h"
RUN_SCRIPT_DIR=${RECOVERY_ROOT}/code/_h
cp "${SOURCE_ROOT}/config/expected-tasks.tsv" "${RECOVERY_ROOT}/config/expected-tasks.tsv"
cp "${SOURCE_ROOT}/config/elastic-net-calibration.rds" \
    "${RECOVERY_ROOT}/config/elastic-net-calibration.rds"
for file in acceptance-criteria.tsv calibration-acceptance-results.tsv \
    calibration-performance-overall.tsv; do
    cp "${SOURCE_ROOT}/config/${file}" "${RECOVERY_ROOT}/config/${file}"
done
CALIBRATION_MODEL=${RECOVERY_ROOT}/config/elastic-net-calibration.rds

job_ids=()
printf 'stage\tregion\tjob_id\tdependency\n' > \
    "${RECOVERY_ROOT}/provenance/submitted-jobs.tsv"
for region in "${REGIONS[@]}"; do
    source_region=${SOURCE_ROOT}/results/${region}/${POPULATION}
    recovery_region=${RECOVERY_ROOT}/results/${region}/${POPULATION}
    mkdir -p "${recovery_region}/summary" "${recovery_region}/excluded" \
        "${recovery_region}/qc_failures"
    cp -a "${source_region}/summary/." "${recovery_region}/summary/"
    for exclusion_file in "${source_region}"/excluded/vmr-*.tsv; do
        [[ -e "${exclusion_file}" ]] || continue
        if ! grep -Fq 'no_variants_in_upstream_plink_window' "${exclusion_file}"; then
            cp "${exclusion_file}" "${recovery_region}/excluded/"
        fi
    done
    if [[ -d "${source_region}/qc_failures" ]]; then
        cp -a "${source_region}/qc_failures/." "${recovery_region}/qc_failures/"
    fi

    manifest=${RECOVERY_ROOT}/config/retry-tasks-${region}.tsv
    printf 'task_id\n' > "${manifest}"
    {
        for failure_file in "${source_region}"/failures/vmr-*.tsv; do
            [[ -e "${failure_file}" ]] || continue
            awk -F '\t' 'FNR > 1 {print $1}' "${failure_file}"
        done
        for exclusion_file in "${source_region}"/excluded/vmr-*.tsv; do
            [[ -e "${exclusion_file}" ]] || continue
            if grep -Fq 'no_variants_in_upstream_plink_window' "${exclusion_file}"; then
                awk -F '\t' 'FNR == 2 {print $1}' "${exclusion_file}"
            fi
        done
    } | sort -n -u >> "${manifest}"
    tasks=$(($(wc -l < "${manifest}") - 1))
    if (( tasks < 1 )); then
        echo "No failed tasks found for ${region}; copied terminal records only"
        continue
    fi
    job=$(sbatch --parsable --account="${ACCOUNT}" \
        --array="1-${tasks}%${MAX_CONCURRENT}" \
        --job-name="cal_h2_recover_${region}" \
        --output="${RECOVERY_ROOT}/logs/%x.%A_%a.log" \
        --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_REPO_ROOT="${REPO_ROOT}",CAL_H2_CALIBRATION_MODEL="${CALIBRATION_MODEL}",CAL_H2_OBSERVED_OUTPUT_ROOT="${RECOVERY_ROOT}/results",CAL_H2_PLINK_ROOT="${PLINK_ROOT}",CAL_H2_RECOVERED_PLINK_ROOT="${RECOVERED_PLINK_ROOT}",CAL_H2_PGEN_ROOT="${PGEN_ROOT}",CAL_H2_PHENOTYPE_ROOT="${PHENOTYPE_ROOT}",CAL_H2_TASK_MANIFEST="${manifest}",CAL_H2_WRITE_DIAGNOSTICS=FALSE,REGION="${region}",POPULATION="${POPULATION}" \
        "${RUN_SCRIPT_DIR}/step_5_estimate_observed_vmr.sh")
    job_id=${job%%;*}
    job_ids+=("${job_id}")
    printf 'recover\t%s\t%s\tNA\n' "${region}" "${job_id}" >> \
        "${RECOVERY_ROOT}/provenance/submitted-jobs.tsv"
done

if (( ${#job_ids[@]} == 0 )); then
    echo "Source run has no computational failures requiring recovery" >&2
    exit 1
fi
dependency=$(IFS=:; echo "${job_ids[*]}")
qc_job=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${dependency}" \
    --output="${RECOVERY_ROOT}/logs/%x.%j.log" \
    --export=ALL,CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_OBSERVED_OUTPUT_ROOT="${RECOVERY_ROOT}/results",CAL_H2_EXPECTED_TASKS="${RECOVERY_ROOT}/config/expected-tasks.tsv" \
    "${RUN_SCRIPT_DIR}/step_6_combine_observed.sh")
qc_job_id=${qc_job%%;*}
printf 'combine_qc\tall\t%s\tafterok:%s\n' "${qc_job_id}" "${dependency}" >> \
    "${RECOVERY_ROOT}/provenance/submitted-jobs.tsv"

{
    printf 'field\tvalue\n'
    printf 'source_observed_run_id\t%s\n' "${SOURCE_RUN_ID}"
    printf 'recovery_run_id\t%s\n' "${RECOVERY_RUN_ID}"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'population\t%s\n' "${POPULATION}"
    printf 'plink_root\t%s\n' "${PLINK_ROOT}"
    printf 'phenotype_fallback_root\t%s\n' "${PHENOTYPE_ROOT}"
    printf 'cohort_pgen_root\t%s\n' "${PGEN_ROOT}"
    printf 'recovered_plink_root\t%s\n' "${RECOVERED_PLINK_ROOT}"
    printf 'recovery_policy\tmissing upstream PLINK inputs are reconstructed from the cohort PGEN; no SNP in +/-500 kb and fewer than two SNPs after QC are QC failures; computational failures have zero tolerance\n'
} > "${RECOVERY_ROOT}/provenance/recovery-metadata.tsv"
git -C "${REPO_ROOT}" rev-parse HEAD > "${RECOVERY_ROOT}/provenance/git-commit.txt"
conda list -p "${ENV_PATH}" --explicit > \
    "${RECOVERY_ROOT}/provenance/conda-explicit-spec.txt"
for path in "${RUN_SCRIPT_DIR}"/*.R "${RUN_SCRIPT_DIR}"/*.sh \
    "${RECOVERY_ROOT}/config"/*; do sha256sum "${path}"; done > \
    "${RECOVERY_ROOT}/provenance/sha256sums.txt"

echo "Recovery run: ${RECOVERY_ROOT}"
echo "Recovery arrays: ${job_ids[*]}"
echo "Combine/QC job: ${qc_job_id}"
