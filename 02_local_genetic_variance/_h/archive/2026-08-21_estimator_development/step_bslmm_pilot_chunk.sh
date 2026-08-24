#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=bslmm_en_pilot
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

SCRIPT_DIR=${CAL_H2_SCRIPT_DIR:?CAL_H2_SCRIPT_DIR must be set}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
MANIFEST=${BSLMM_PILOT_MANIFEST:?BSLMM_PILOT_MANIFEST must be set}
CHUNK_MANIFEST=${BSLMM_PILOT_CHUNK_MANIFEST:?BSLMM_PILOT_CHUNK_MANIFEST must be set}
OUTPUT_ROOT=${BSLMM_PILOT_OUTPUT_ROOT:?BSLMM_PILOT_OUTPUT_ROOT must be set}
CALIBRATION_MODEL=${BSLMM_PILOT_CALIBRATION_MODEL:?}
EXPECTED_SHA=${BSLMM_PILOT_CALIBRATION_SHA256:?}
WORK_ROOT=${BSLMM_PILOT_WORK_ROOT:-${OUTPUT_ROOT}/work}
KEEP_WORK=${BSLMM_PILOT_KEEP_WORK:-FALSE}

export PATH="/projects/p32505/opt/bin:${PATH}"

"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/13_run_bslmm_en_pilot_scenario.R" \
    --manifest="${MANIFEST}" \
    --chunk_manifest="${CHUNK_MANIFEST}" \
    --chunk_id="${SLURM_ARRAY_TASK_ID}" \
    --calibration_model="${CALIBRATION_MODEL}" \
    --expected_calibration_sha256="${EXPECTED_SHA}" \
    --output_root="${OUTPUT_ROOT}/raw" \
    --work_root="${WORK_ROOT}" \
    --keep_work="${KEEP_WORK}"
