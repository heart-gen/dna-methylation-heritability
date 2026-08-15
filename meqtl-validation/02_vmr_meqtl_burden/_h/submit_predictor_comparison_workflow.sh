#!/bin/bash

set -euo pipefail

usage() {
    echo "Usage: $0 RUN_ID"
}
if (( $# != 1 )); then
    usage >&2
    exit 1
fi

RUN_ID=$1
if [[ ! "${RUN_ID}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "RUN_ID contains unsupported characters: ${RUN_ID}" >&2
    exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MODULE_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${MODULE_DIR}/../.." && pwd)
RUN_ROOT=${MODULE_DIR}/_m/predictor-comparison-runs/${RUN_ID}
ENV_PATH=${PREDICTOR_COMPARISON_ENV:-/projects/p32505/opt/envs/genomics}
ACCOUNT=${SBATCH_ACCOUNT:-b1042}
PARTITION=${SBATCH_PARTITION:-genomics}
CONFIG_SOURCE=${REPO_ROOT}/config/predictor_comparison.yml
CALIBRATED_ROOT=${REPO_ROOT}/calibrated-simulation-analysis/_m/observed-runs/observed-AA-v4/results/combined
CALIBRATED=${CALIBRATED_ROOT}/calibrated-local-h2-AA-vmrs.tsv
CALIBRATED_QC=${CALIBRATED_ROOT}/observed-run-qc.tsv
REGIONS=(caudate dlpfc hippocampus)

if [[ -e "${RUN_ROOT}" ]]; then
    echo "Run already exists and will not be overwritten: ${RUN_ROOT}" >&2
    exit 1
fi
for path in "${CONFIG_SOURCE}" "${CALIBRATED}" "${CALIBRATED_QC}" \
    "${ENV_PATH}/bin/python"; do
    if [[ ! -s "${path}" ]]; then
        echo "Required input is missing or empty: ${path}" >&2
        exit 1
    fi
done
for region in "${REGIONS[@]}"; do
    for path in \
        "${REPO_ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${region}/_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz" \
        "${REPO_ROOT}/heritability/elastic_net_model/all_individuals/${region}/_m/${region}_summary_elastic-net_AA.tsv" \
        "${REPO_ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_m/${region}/vmr_technical_annotations.tsv"; do
        if [[ ! -s "${path}" ]]; then
            echo "Required region input is missing or empty: ${path}" >&2
            exit 1
        fi
    done
done

"${ENV_PATH}/bin/python" -c \
    'import matplotlib, numpy, pandas, scipy, statsmodels, yaml; print("comparison environment verified")'
"${ENV_PATH}/bin/python" - "${CALIBRATED_QC}" <<'PY'
import sys
import pandas as pd

path = sys.argv[1]
qc = pd.read_csv(path, sep="\t")
expected = {"caudate", "dlpfc", "hippocampus"}
observed = set(qc["region"].astype(str).str.lower())
truth = qc["overall_qc_pass"].astype(str).str.lower().isin({"true", "1", "yes"})
if observed != expected or not truth.all() or not qc["computational_failed_tasks"].eq(0).all():
    raise SystemExit("Accepted observed-AA-v4 QC gate did not pass")
print("accepted calibrated run verified")
PY

mkdir -p "${RUN_ROOT}/code/meqtl-validation/02_vmr_meqtl_burden" \
    "${RUN_ROOT}/code/meqtl-validation/07_repeat_mappability_sensitivity/_h" \
    "${RUN_ROOT}/code/meqtl-validation/_lib" "${RUN_ROOT}/code/config" \
    "${RUN_ROOT}/config" "${RUN_ROOT}/logs" "${RUN_ROOT}/provenance"
cp -a "${SCRIPT_DIR}" "${RUN_ROOT}/code/meqtl-validation/02_vmr_meqtl_burden/_h"
cp -a "${REPO_ROOT}/meqtl-validation/_lib/." \
    "${RUN_ROOT}/code/meqtl-validation/_lib/"
cp "${REPO_ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_h/04_complete_tech_joins.py" \
    "${RUN_ROOT}/code/meqtl-validation/07_repeat_mappability_sensitivity/_h/"
cp "${REPO_ROOT}/config/paths.yml" "${REPO_ROOT}/config/meqtl_parameters.yml" \
    "${REPO_ROOT}/config/analysis_thresholds.yml" "${CONFIG_SOURCE}" \
    "${RUN_ROOT}/code/config/"
cp "${CONFIG_SOURCE}" "${RUN_ROOT}/config/predictor_comparison.yml"
cp "${CALIBRATED}" "${RUN_ROOT}/config/calibrated-local-h2-AA-vmrs.tsv"
cp "${CALIBRATED_QC}" "${RUN_ROOT}/config/observed-run-qc.tsv"
CODE_ROOT=${RUN_ROOT}/code
RUN_SCRIPT_DIR=${CODE_ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h

{
    printf 'field\tvalue\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'population\tAA\n'
    printf 'calibrated_source_run\tobserved-AA-v4\n'
    printf 'environment\t%s\n' "${ENV_PATH}"
    printf 'account\t%s\n' "${ACCOUNT}"
    printf 'partition\t%s\n' "${PARTITION}"
} > "${RUN_ROOT}/provenance/run-metadata.tsv"
git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_ROOT}/provenance/git-commit.txt"
conda list -p "${ENV_PATH}" --explicit > "${RUN_ROOT}/provenance/conda-explicit-spec.txt"

INPUT_MANIFEST=${RUN_ROOT}/provenance/input-manifest.tsv
printf 'input_type\tregion\tpath\tbytes\tsha256\n' > "${INPUT_MANIFEST}"
record_input() {
    local input_type=$1
    local region=$2
    local path=$3
    printf '%s\t%s\t%s\t%s\t%s\n' "${input_type}" "${region}" "${path}" \
        "$(stat -c %s "${path}")" "$(sha256sum "${path}" | cut -d ' ' -f 1)" \
        >> "${INPUT_MANIFEST}"
}
record_input calibrated_estimates all "${CALIBRATED}"
record_input calibrated_qc all "${CALIBRATED_QC}"
record_input comparison_config all "${CONFIG_SOURCE}"
for region in "${REGIONS[@]}"; do
    record_input phase1_m3a_leads "${region}" \
        "${REPO_ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${region}/_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz"
    record_input legacy_predictability "${region}" \
        "${REPO_ROOT}/heritability/elastic_net_model/all_individuals/${region}/_m/${region}_summary_elastic-net_AA.tsv"
    record_input technical_annotations "${region}" \
        "${REPO_ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_m/${region}/vmr_technical_annotations.tsv"
done
find "${RUN_ROOT}/code" "${RUN_ROOT}/config" -type f -print0 | sort -z | \
    xargs -0 sha256sum > "${RUN_ROOT}/provenance/code-config-sha256sums.txt"

export_args=ALL,PREDICTOR_COMPARISON_RUN_ROOT="${RUN_ROOT}",PREDICTOR_COMPARISON_CODE_ROOT="${CODE_ROOT}",PREDICTOR_COMPARISON_REPO_ROOT="${REPO_ROOT}",PREDICTOR_COMPARISON_ENV="${ENV_PATH}"
printf 'stage\tjob_id\tdependency\n' > "${RUN_ROOT}/provenance/submitted-jobs.tsv"
aggregate=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --output="${RUN_ROOT}/logs/%x.%A_%a.log" --export="${export_args}" \
    "${RUN_SCRIPT_DIR}/step_3a_aggregate_repaired.sh")
aggregate_id=${aggregate%%;*}
printf 'aggregate\t%s\tNA\n' "${aggregate_id}" >> "${RUN_ROOT}/provenance/submitted-jobs.tsv"

tech_join=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --dependency="afterok:${aggregate_id}" --output="${RUN_ROOT}/logs/%x.%j.log" \
    --export="${export_args}" "${RUN_SCRIPT_DIR}/step_3b_repair_technical_join.sh")
tech_join_id=${tech_join%%;*}
printf 'technical_join\t%s\tafterok:%s\n' "${tech_join_id}" "${aggregate_id}" \
    >> "${RUN_ROOT}/provenance/submitted-jobs.tsv"

compare=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --dependency="afterok:${tech_join_id}" --output="${RUN_ROOT}/logs/%x.%A_%a.log" \
    --export="${export_args}" "${RUN_SCRIPT_DIR}/step_3c_compare_predictors.sh")
compare_id=${compare%%;*}
printf 'comparison\t%s\tafterok:%s\n' "${compare_id}" "${tech_join_id}" \
    >> "${RUN_ROOT}/provenance/submitted-jobs.tsv"

combine=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --dependency="afterok:${compare_id}" --output="${RUN_ROOT}/logs/%x.%j.log" \
    --export="${export_args}" "${RUN_SCRIPT_DIR}/step_3d_combine_and_plot.sh")
combine_id=${combine%%;*}
printf 'combine_plot\t%s\tafterok:%s\n' "${combine_id}" "${compare_id}" \
    >> "${RUN_ROOT}/provenance/submitted-jobs.tsv"

echo "Predictor-comparison run: ${RUN_ROOT}"
echo "Jobs: aggregate=${aggregate_id}, technical_join=${tech_join_id}, comparison=${compare_id}, combine_plot=${combine_id}"

