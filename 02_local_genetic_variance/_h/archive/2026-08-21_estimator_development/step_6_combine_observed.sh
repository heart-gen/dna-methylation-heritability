#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_obs_qc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00

set -euo pipefail

SCRIPT_DIR=${CAL_H2_SCRIPT_DIR:?CAL_H2_SCRIPT_DIR must be set}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
OUTPUT_ROOT=${CAL_H2_OBSERVED_OUTPUT_ROOT:?CAL_H2_OBSERVED_OUTPUT_ROOT must be set}
EXPECTED=${CAL_H2_EXPECTED_TASKS:?CAL_H2_EXPECTED_TASKS must be set}

"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/07_combine_observed.R" \
    --input="${OUTPUT_ROOT}" \
    --expected="${EXPECTED}" \
    --output-dir="${OUTPUT_ROOT}/combined" \
    --fail-on-qc=TRUE
