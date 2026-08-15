#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=predfu_extfix

set -euo pipefail

: "${PREDICTOR_FOLLOWUP_RUN_ROOT:?PREDICTOR_FOLLOWUP_RUN_ROOT must be set}"
: "${PREDICTOR_FOLLOWUP_CODE_ROOT:?PREDICTOR_FOLLOWUP_CODE_ROOT must be set}"
: "${PREDICTOR_FOLLOWUP_REPO_ROOT:?PREDICTOR_FOLLOWUP_REPO_ROOT must be set}"

ENV_PATH=${PREDICTOR_FOLLOWUP_ENV:-/projects/p32505/opt/envs/genomics}
SCRIPT_DIR=${PREDICTOR_FOLLOWUP_CODE_ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h
REGION=dlpfc

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/11_fit_orthogonal_validation.py" \
    --region "${REGION}" \
    --primary-bridge "${PREDICTOR_FOLLOWUP_RUN_ROOT}/inputs/base-run/${REGION}/predictor_bridge.tsv.gz" \
    --calibrated "${PREDICTOR_FOLLOWUP_RUN_ROOT}/config/calibrated-local-h2-AA-vmrs.tsv" \
    --external-root "${PREDICTOR_FOLLOWUP_REPO_ROOT}/meqtl-validation/03_external_meqtl_validation/_m" \
    --regulatory-context-root "${PREDICTOR_FOLLOWUP_REPO_ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/regulatory_context/_m" \
    --output-dir "${PREDICTOR_FOLLOWUP_RUN_ROOT}/orthogonal/${REGION}"

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/12_combine_followup.py" \
    --input-root "${PREDICTOR_FOLLOWUP_RUN_ROOT}" \
    --followup-config "${PREDICTOR_FOLLOWUP_RUN_ROOT}/config/predictor_followup.yml" \
    --base-promotion "${PREDICTOR_FOLLOWUP_RUN_ROOT}/inputs/base-run/promotion_status.tsv" \
    --output-dir "${PREDICTOR_FOLLOWUP_RUN_ROOT}/combined"

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/13_plot_followup.py" \
    --base-combined "${PREDICTOR_FOLLOWUP_RUN_ROOT}/inputs/base-run/combined" \
    --followup-combined "${PREDICTOR_FOLLOWUP_RUN_ROOT}/combined" \
    --output-dir "${PREDICTOR_FOLLOWUP_RUN_ROOT}/figures"

