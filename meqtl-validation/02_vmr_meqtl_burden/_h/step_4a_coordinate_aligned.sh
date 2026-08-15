#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --array=0-2
#SBATCH --job-name=predfu_align

set -euo pipefail

: "${PREDICTOR_FOLLOWUP_RUN_ROOT:?PREDICTOR_FOLLOWUP_RUN_ROOT must be set}"
: "${PREDICTOR_FOLLOWUP_CODE_ROOT:?PREDICTOR_FOLLOWUP_CODE_ROOT must be set}"
: "${PREDICTOR_FOLLOWUP_REPO_ROOT:?PREDICTOR_FOLLOWUP_REPO_ROOT must be set}"

REGIONS=(caudate dlpfc hippocampus)
REGION=${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}
ENV_PATH=${PREDICTOR_FOLLOWUP_ENV:-/projects/p32505/opt/envs/genomics}
SCRIPT_DIR=${PREDICTOR_FOLLOWUP_CODE_ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h
OUTDIR=${PREDICTOR_FOLLOWUP_RUN_ROOT}/coordinate_aligned/${REGION}
CALIBRATED=${PREDICTOR_FOLLOWUP_RUN_ROOT}/config/calibrated-local-h2-AA-vmrs.tsv
CALIBRATED_QC=${PREDICTOR_FOLLOWUP_RUN_ROOT}/config/observed-run-qc.tsv
PREDICTOR_CONFIG=${PREDICTOR_FOLLOWUP_RUN_ROOT}/config/predictor_comparison.yml
LEGACY=${PREDICTOR_FOLLOWUP_REPO_ROOT}/heritability/elastic_net_model/all_individuals/${REGION}/_m/${REGION}_summary_elastic-net_AA.tsv
TECHNICAL=${PREDICTOR_FOLLOWUP_REPO_ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_m/${REGION}/vmr_technical_annotations.tsv

mkdir -p "${OUTDIR}"
"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/08_reaggregate_coordinate_aligned.py" \
    --region "${REGION}" \
    --lead "${PREDICTOR_FOLLOWUP_REPO_ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz" \
    --vmr-bed "${PREDICTOR_FOLLOWUP_REPO_ROOT}/vmr-analysis/all_individuals/${REGION}/_m/vmr.bed" \
    --prepared-map-glob "${PREDICTOR_FOLLOWUP_REPO_ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/prepared/cpg_vmr_map.chr*.tsv" \
    --legacy "${LEGACY}" \
    --technical "${TECHNICAL}" \
    --annotation "${PREDICTOR_FOLLOWUP_REPO_ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/annotation/_m/${REGION}_vmr_annotations_hg38.tsv" \
    --output-dir "${OUTDIR}" \
    --fdr 0.05

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/04_build_predictor_bridge.py" \
    --region "${REGION}" \
    --burden "${OUTDIR}/vmr_meqtl_burden.tsv.gz" \
    --legacy "${LEGACY}" \
    --calibrated "${CALIBRATED}" \
    --calibrated-qc "${CALIBRATED_QC}" \
    --config "${PREDICTOR_CONFIG}" \
    --output-dir "${OUTDIR}"

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/05_fit_predictor_comparison.py" \
    --region "${REGION}" \
    --bridge "${OUTDIR}/predictor_bridge.tsv.gz" \
    --config "${PREDICTOR_CONFIG}" \
    --output-dir "${OUTDIR}"

