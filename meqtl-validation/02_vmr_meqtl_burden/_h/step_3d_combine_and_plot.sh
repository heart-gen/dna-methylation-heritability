#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=predcmp_combine

set -euo pipefail

: "${PREDICTOR_COMPARISON_RUN_ROOT:?PREDICTOR_COMPARISON_RUN_ROOT must be set}"
: "${PREDICTOR_COMPARISON_CODE_ROOT:?PREDICTOR_COMPARISON_CODE_ROOT must be set}"

ENV_PATH=${PREDICTOR_COMPARISON_ENV:-/projects/p32505/opt/envs/genomics}
SCRIPT_DIR=${PREDICTOR_COMPARISON_CODE_ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h
CONFIG=${PREDICTOR_COMPARISON_RUN_ROOT}/config/predictor_comparison.yml
COMBINED=${PREDICTOR_COMPARISON_RUN_ROOT}/combined
FIGURES=${PREDICTOR_COMPARISON_RUN_ROOT}/figures

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/06_combine_predictor_comparison.py" \
    --input-root "${PREDICTOR_COMPARISON_RUN_ROOT}/regions" \
    --config "${CONFIG}" \
    --output-dir "${COMBINED}"

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/07_plot_predictor_comparison.py" \
    --region-root "${PREDICTOR_COMPARISON_RUN_ROOT}/regions" \
    --combined-dir "${COMBINED}" \
    --output-dir "${FIGURES}" \
    --seed 20260808

