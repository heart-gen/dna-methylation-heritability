#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --array=0-2
#SBATCH --job-name=predfu_sens

set -euo pipefail

: "${PREDICTOR_FOLLOWUP_RUN_ROOT:?PREDICTOR_FOLLOWUP_RUN_ROOT must be set}"
: "${PREDICTOR_FOLLOWUP_CODE_ROOT:?PREDICTOR_FOLLOWUP_CODE_ROOT must be set}"
: "${PREDICTOR_FOLLOWUP_REPO_ROOT:?PREDICTOR_FOLLOWUP_REPO_ROOT must be set}"

REGIONS=(caudate dlpfc hippocampus)
REGION=${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}
ENV_PATH=${PREDICTOR_FOLLOWUP_ENV:-/projects/p32505/opt/envs/genomics}
SCRIPT_DIR=${PREDICTOR_FOLLOWUP_CODE_ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h
PREDICTOR_CONFIG=${PREDICTOR_FOLLOWUP_RUN_ROOT}/config/predictor_comparison.yml
CALIBRATED=${PREDICTOR_FOLLOWUP_RUN_ROOT}/config/calibrated-local-h2-AA-vmrs.tsv
PRIMARY_BRIDGE=${PREDICTOR_FOLLOWUP_RUN_ROOT}/inputs/base-run/${REGION}/predictor_bridge.tsv.gz
ALIGNED_BRIDGE=${PREDICTOR_FOLLOWUP_RUN_ROOT}/coordinate_aligned/${REGION}/predictor_bridge.tsv.gz
ANNOTATION_ROOT=${PREDICTOR_FOLLOWUP_REPO_ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/annotation
TECHNICAL=${PREDICTOR_FOLLOWUP_REPO_ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_m/${REGION}/vmr_technical_annotations.tsv

mkdir -p "${PREDICTOR_FOLLOWUP_RUN_ROOT}/sensitivities/${REGION}" \
    "${PREDICTOR_FOLLOWUP_RUN_ROOT}/annotation/${REGION}" \
    "${PREDICTOR_FOLLOWUP_RUN_ROOT}/orthogonal/${REGION}"

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/09_fit_boundary_and_weighting.py" \
    --region "${REGION}" \
    --primary-bridge "${PRIMARY_BRIDGE}" \
    --aligned-bridge "${ALIGNED_BRIDGE}" \
    --config "${PREDICTOR_CONFIG}" \
    --output-dir "${PREDICTOR_FOLLOWUP_RUN_ROOT}/sensitivities/${REGION}"

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/10_fit_annotation_sensitivity.py" \
    --region "${REGION}" \
    --repeat "${ANNOTATION_ROOT}/repeat_elements/_m/vmr_repeat_overlap_AA.tsv" \
    --repressive "${ANNOTATION_ROOT}/repressive_chromatin/_m/vmr_repressive_overlap_AA.tsv" \
    --technical "${TECHNICAL}" \
    --calibrated "${CALIBRATED}" \
    --output-dir "${PREDICTOR_FOLLOWUP_RUN_ROOT}/annotation/${REGION}" \
    --seed 20260809 \
    --n-permutations 10000

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/11_fit_orthogonal_validation.py" \
    --region "${REGION}" \
    --primary-bridge "${PRIMARY_BRIDGE}" \
    --calibrated "${CALIBRATED}" \
    --external-root "${PREDICTOR_FOLLOWUP_REPO_ROOT}/meqtl-validation/03_external_meqtl_validation/_m" \
    --regulatory-context-root "${PREDICTOR_FOLLOWUP_REPO_ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/regulatory_context/_m" \
    --output-dir "${PREDICTOR_FOLLOWUP_RUN_ROOT}/orthogonal/${REGION}"

