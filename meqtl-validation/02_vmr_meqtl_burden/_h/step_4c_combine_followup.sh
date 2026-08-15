#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=predfu_combine

set -euo pipefail

: "${PREDICTOR_FOLLOWUP_RUN_ROOT:?PREDICTOR_FOLLOWUP_RUN_ROOT must be set}"
: "${PREDICTOR_FOLLOWUP_CODE_ROOT:?PREDICTOR_FOLLOWUP_CODE_ROOT must be set}"

ENV_PATH=${PREDICTOR_FOLLOWUP_ENV:-/projects/p32505/opt/envs/genomics}
SCRIPT_DIR=${PREDICTOR_FOLLOWUP_CODE_ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h
COMBINED=${PREDICTOR_FOLLOWUP_RUN_ROOT}/combined

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/12_combine_followup.py" \
    --input-root "${PREDICTOR_FOLLOWUP_RUN_ROOT}" \
    --followup-config "${PREDICTOR_FOLLOWUP_RUN_ROOT}/config/predictor_followup.yml" \
    --base-promotion "${PREDICTOR_FOLLOWUP_RUN_ROOT}/inputs/base-run/promotion_status.tsv" \
    --output-dir "${COMBINED}"

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/13_plot_followup.py" \
    --base-combined "${PREDICTOR_FOLLOWUP_RUN_ROOT}/inputs/base-run/combined" \
    --followup-combined "${COMBINED}" \
    --output-dir "${PREDICTOR_FOLLOWUP_RUN_ROOT}/figures"

