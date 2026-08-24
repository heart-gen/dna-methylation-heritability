#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=joint_pve
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

SCRIPT_DIR=${JOINT_PVE_SCRIPT_DIR:?}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
MANIFEST=${JOINT_PVE_MANIFEST:?}
CHUNK_MANIFEST=${JOINT_PVE_CHUNK_MANIFEST:?}
OUTPUT_ROOT=${JOINT_PVE_OUTPUT_ROOT:?}
WORK_ROOT=${JOINT_PVE_WORK_ROOT:-${OUTPUT_ROOT}/work}
DEVELOPMENT_BSLMM_ROOT=${JOINT_PVE_DEVELOPMENT_BSLMM_ROOT:-}
KEEP_WORK=${JOINT_PVE_KEEP_WORK:-FALSE}

export PATH="/projects/p32505/opt/bin:${PATH}"

"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/22_run_joint_pve_scenario.R" \
    --manifest="${MANIFEST}" \
    --chunk_manifest="${CHUNK_MANIFEST}" \
    --chunk_id="${SLURM_ARRAY_TASK_ID}" \
    --output_root="${OUTPUT_ROOT}/raw" \
    --work_root="${WORK_ROOT}" \
    --development_bslmm_root="${DEVELOPMENT_BSLMM_ROOT}" \
    --keep_work="${KEEP_WORK}"
