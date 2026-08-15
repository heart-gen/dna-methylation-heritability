#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --job-name=predcmp_techjoin

set -euo pipefail

: "${PREDICTOR_COMPARISON_RUN_ROOT:?PREDICTOR_COMPARISON_RUN_ROOT must be set}"
: "${PREDICTOR_COMPARISON_CODE_ROOT:?PREDICTOR_COMPARISON_CODE_ROOT must be set}"
: "${PREDICTOR_COMPARISON_REPO_ROOT:?PREDICTOR_COMPARISON_REPO_ROOT must be set}"

ENV_PATH=${PREDICTOR_COMPARISON_ENV:-/projects/p32505/opt/envs/genomics}
SCRIPT=${PREDICTOR_COMPARISON_CODE_ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_h/04_complete_tech_joins.py

"${ENV_PATH}/bin/python" "${SCRIPT}" \
    --join-only \
    --phase2-root "${PREDICTOR_COMPARISON_RUN_ROOT}/regions" \
    --technical-root "${PREDICTOR_COMPARISON_REPO_ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_m" \
    --join-output-root "${PREDICTOR_COMPARISON_RUN_ROOT}/technical-joins" \
    --elastic-net-root "${PREDICTOR_COMPARISON_REPO_ROOT}/heritability/elastic_net_model/all_individuals" \
    --min-reciprocal-overlap 0.5

