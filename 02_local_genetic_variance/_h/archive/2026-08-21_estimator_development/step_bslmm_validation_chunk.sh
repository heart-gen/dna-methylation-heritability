#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=bslmm_val
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

SCRIPT_DIR=${CAL_H2_SCRIPT_DIR:?}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
MANIFEST=${BSLMM_VAL_MANIFEST:?}
CHUNK_MANIFEST=${BSLMM_VAL_CHUNK_MANIFEST:?}
OUTPUT_ROOT=${BSLMM_VAL_OUTPUT_ROOT:?}
WORK_ROOT=${BSLMM_VAL_WORK_ROOT:-${OUTPUT_ROOT}/work}
KEEP_WORK=${BSLMM_VAL_KEEP_WORK:-FALSE}

export PATH="/projects/p32505/opt/bin:${PATH}"

"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/16_run_bslmm_validation_scenario.R" \
    --manifest="${MANIFEST}" \
    --chunk_manifest="${CHUNK_MANIFEST}" \
    --chunk_id="${SLURM_ARRAY_TASK_ID}" \
    --output_root="${OUTPUT_ROOT}/raw" \
    --work_root="${WORK_ROOT}" \
    --keep_work="${KEEP_WORK}"
