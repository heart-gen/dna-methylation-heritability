#!/bin/bash

set -euo pipefail

if (( $# != 1 )); then
    echo "Usage: $0 RUN_ID" >&2
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
RUN_ROOT=${MODULE_DIR}/_m/predictor-followup-runs/${RUN_ID}
BASE_ROOT=${MODULE_DIR}/_m/predictor-comparison-runs/calibrated-vs-legacy-AA-v1
ENV_PATH=${PREDICTOR_FOLLOWUP_ENV:-/projects/p32505/opt/envs/genomics}
ACCOUNT=${SBATCH_ACCOUNT:-b1042}
PARTITION=${SBATCH_PARTITION:-genomics}
CALIBRATED_ROOT=${REPO_ROOT}/calibrated-simulation-analysis/_m/observed-runs/observed-AA-v4/results/combined
CALIBRATED=${CALIBRATED_ROOT}/calibrated-local-h2-AA-vmrs.tsv
CALIBRATED_QC=${CALIBRATED_ROOT}/observed-run-qc.tsv
REGIONS=(caudate dlpfc hippocampus)

if [[ -e "${RUN_ROOT}" ]]; then
    echo "Run already exists and will not be overwritten: ${RUN_ROOT}" >&2
    exit 1
fi
for path in \
    "${ENV_PATH}/bin/python" \
    "${CALIBRATED}" "${CALIBRATED_QC}" \
    "${BASE_ROOT}/combined/promotion_status.tsv" \
    "${REPO_ROOT}/config/predictor_comparison.yml" \
    "${REPO_ROOT}/config/predictor_followup.yml"; do
    if [[ ! -s "${path}" ]]; then
        echo "Required input is missing or empty: ${path}" >&2
        exit 1
    fi
done
for region in "${REGIONS[@]}"; do
    for path in \
        "${BASE_ROOT}/regions/${region}/predictor_bridge.tsv.gz" \
        "${REPO_ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${region}/_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz" \
        "${REPO_ROOT}/vmr-analysis/all_individuals/${region}/_m/vmr.bed" \
        "${REPO_ROOT}/heritability/elastic_net_model/all_individuals/${region}/_m/${region}_summary_elastic-net_AA.tsv" \
        "${REPO_ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_m/${region}/vmr_technical_annotations.tsv"; do
        if [[ ! -s "${path}" ]]; then
            echo "Required region input is missing or empty: ${path}" >&2
            exit 1
        fi
    done
done

"${ENV_PATH}/bin/python" -c \
    'import matplotlib, numpy, pandas, scipy, statsmodels, yaml; print("follow-up environment verified")'
"${ENV_PATH}/bin/python" - "${BASE_ROOT}/combined/promotion_status.tsv" "${CALIBRATED_QC}" <<'PY'
import sys
import pandas as pd

promotion = pd.read_csv(sys.argv[1], sep="\t")
qc = pd.read_csv(sys.argv[2], sep="\t")
as_bool = lambda values: values.astype(str).str.lower().isin({"true", "1", "yes"})
if len(promotion) != 1 or not as_bool(promotion["phase2_calibrated_gate_pass"]).iloc[0]:
    raise SystemExit("Base Phase 2 calibrated gate is not accepted")
if set(qc["region"].str.lower()) != {"caudate", "dlpfc", "hippocampus"}:
    raise SystemExit("Calibrated QC does not contain exactly the three expected regions")
if not as_bool(qc["overall_qc_pass"]).all() or not qc["computational_failed_tasks"].eq(0).all():
    raise SystemExit("Accepted calibrated observed-run QC failed")
print("base and calibrated acceptance gates verified")
PY

mkdir -p "${RUN_ROOT}/code/meqtl-validation/02_vmr_meqtl_burden" \
    "${RUN_ROOT}/code/meqtl-validation/_lib" "${RUN_ROOT}/code/config" \
    "${RUN_ROOT}/config" "${RUN_ROOT}/inputs/base-run/combined" \
    "${RUN_ROOT}/logs" "${RUN_ROOT}/provenance"
cp -a "${SCRIPT_DIR}" "${RUN_ROOT}/code/meqtl-validation/02_vmr_meqtl_burden/_h"
cp -a "${REPO_ROOT}/meqtl-validation/_lib/." "${RUN_ROOT}/code/meqtl-validation/_lib/"
cp "${REPO_ROOT}/config/paths.yml" "${REPO_ROOT}/config/meqtl_parameters.yml" \
    "${REPO_ROOT}/config/analysis_thresholds.yml" \
    "${REPO_ROOT}/config/predictor_comparison.yml" \
    "${REPO_ROOT}/config/predictor_followup.yml" "${RUN_ROOT}/code/config/"
cp "${REPO_ROOT}/config/predictor_comparison.yml" "${RUN_ROOT}/config/"
cp "${REPO_ROOT}/config/predictor_followup.yml" "${RUN_ROOT}/config/"
cp "${CALIBRATED}" "${RUN_ROOT}/config/calibrated-local-h2-AA-vmrs.tsv"
cp "${CALIBRATED_QC}" "${RUN_ROOT}/config/observed-run-qc.tsv"
cp -a "${BASE_ROOT}/combined/." "${RUN_ROOT}/inputs/base-run/combined/"
cp "${BASE_ROOT}/combined/promotion_status.tsv" "${RUN_ROOT}/inputs/base-run/promotion_status.tsv"
for region in "${REGIONS[@]}"; do
    mkdir -p "${RUN_ROOT}/inputs/base-run/${region}"
    cp "${BASE_ROOT}/regions/${region}/predictor_bridge.tsv.gz" \
        "${RUN_ROOT}/inputs/base-run/${region}/"
done

CODE_ROOT=${RUN_ROOT}/code
RUN_SCRIPT_DIR=${CODE_ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h
{
    printf 'field\tvalue\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'base_run_id\tcalibrated-vs-legacy-AA-v1\n'
    printf 'calibrated_source_run\tobserved-AA-v4\n'
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'population\tAA\n'
    printf 'environment\t%s\n' "${ENV_PATH}"
    printf 'account\t%s\n' "${ACCOUNT}"
    printf 'partition\t%s\n' "${PARTITION}"
} > "${RUN_ROOT}/provenance/run-metadata.tsv"
git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_ROOT}/provenance/git-commit.txt"
conda list -p "${ENV_PATH}" --explicit > "${RUN_ROOT}/provenance/conda-explicit-spec.txt"

INPUT_MANIFEST=${RUN_ROOT}/provenance/input-manifest.tsv
printf 'input_type\tregion\tpath\tbytes\tsha256\n' > "${INPUT_MANIFEST}"
record_input() {
    local input_type=$1 region=$2 path=$3
    printf '%s\t%s\t%s\t%s\t%s\n' "${input_type}" "${region}" "${path}" \
        "$(stat -c %s "${path}")" "$(sha256sum "${path}" | cut -d ' ' -f 1)" \
        >> "${INPUT_MANIFEST}"
}
record_input calibrated_estimates all "${CALIBRATED}"
record_input calibrated_qc all "${CALIBRATED_QC}"
record_input base_promotion all "${BASE_ROOT}/combined/promotion_status.tsv"
record_input repeat_annotation all "${REPO_ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/annotation/repeat_elements/_m/vmr_repeat_overlap_AA.tsv"
record_input repressive_annotation all "${REPO_ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/annotation/repressive_chromatin/_m/vmr_repressive_overlap_AA.tsv"
record_input jaffe_support dlpfc "${REPO_ROOT}/meqtl-validation/03_external_meqtl_validation/_m/harmonized/jaffe_dlpfc_450k_meqtl.dlpfc.vmr_support.tsv.gz"
record_input schulz_support hippocampus "${REPO_ROOT}/meqtl-validation/03_external_meqtl_validation/_m/harmonized/schulz_hippocampus_array_meqtl.hippocampus.vmr_support.tsv.gz"
for region in "${REGIONS[@]}"; do
    record_input base_bridge "${region}" "${BASE_ROOT}/regions/${region}/predictor_bridge.tsv.gz"
    record_input m3a_leads "${region}" "${REPO_ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${region}/_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz"
    record_input all_individual_vmr "${region}" "${REPO_ROOT}/vmr-analysis/all_individuals/${region}/_m/vmr.bed"
    record_input legacy_predictability "${region}" "${REPO_ROOT}/heritability/elastic_net_model/all_individuals/${region}/_m/${region}_summary_elastic-net_AA.tsv"
    record_input technical_annotations "${region}" "${REPO_ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_m/${region}/vmr_technical_annotations.tsv"
    for modality_path in \
        "${REPO_ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/regulatory_context/_m/${region}/AA/expression/nearest_gene_window_250kb/architecture_model_input.tsv" \
        "${REPO_ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/regulatory_context/_m/${region}/AA/psi/window_250kb/architecture_model_input.tsv" \
        "${REPO_ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/regulatory_context/_m/${region}/AA/expression/abc/architecture_model_input.tsv"; do
        if [[ -s "${modality_path}" ]]; then record_input transcription_input "${region}" "${modality_path}"; fi
    done
done
find "${RUN_ROOT}/code" "${RUN_ROOT}/config" "${RUN_ROOT}/inputs/base-run" -type f -print0 | \
    sort -z | xargs -0 sha256sum > "${RUN_ROOT}/provenance/snapshot-sha256sums.txt"

export_args=ALL,PREDICTOR_FOLLOWUP_RUN_ROOT="${RUN_ROOT}",PREDICTOR_FOLLOWUP_CODE_ROOT="${CODE_ROOT}",PREDICTOR_FOLLOWUP_REPO_ROOT="${REPO_ROOT}",PREDICTOR_FOLLOWUP_ENV="${ENV_PATH}"
printf 'stage\tjob_id\tdependency\n' > "${RUN_ROOT}/provenance/submitted-jobs.tsv"
align=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --output="${RUN_ROOT}/logs/%x.%A_%a.log" --export="${export_args}" \
    "${RUN_SCRIPT_DIR}/step_4a_coordinate_aligned.sh")
align_id=${align%%;*}
printf 'coordinate_aligned\t%s\tNA\n' "${align_id}" >> "${RUN_ROOT}/provenance/submitted-jobs.tsv"

sens=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --dependency="afterok:${align_id}" --output="${RUN_ROOT}/logs/%x.%A_%a.log" \
    --export="${export_args}" "${RUN_SCRIPT_DIR}/step_4b_followup_sensitivities.sh")
sens_id=${sens%%;*}
printf 'sensitivities\t%s\tafterok:%s\n' "${sens_id}" "${align_id}" \
    >> "${RUN_ROOT}/provenance/submitted-jobs.tsv"

combine=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --dependency="afterok:${sens_id}" --output="${RUN_ROOT}/logs/%x.%j.log" \
    --export="${export_args}" "${RUN_SCRIPT_DIR}/step_4c_combine_followup.sh")
combine_id=${combine%%;*}
printf 'combine_plot\t%s\tafterok:%s\n' "${combine_id}" "${sens_id}" \
    >> "${RUN_ROOT}/provenance/submitted-jobs.tsv"

echo "Predictor follow-up run: ${RUN_ROOT}"
echo "Jobs: coordinate=${align_id}, sensitivities=${sens_id}, combine_plot=${combine_id}"

