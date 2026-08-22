#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=lgv_features
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

RUN_DIR=${LGV_RUN_DIR:?LGV_RUN_DIR is required}
H_DIR=${LGV_H_DIR:?LGV_H_DIR is required}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
CHUNK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
CHUNK_MANIFEST=${RUN_DIR}/config/chunk-manifest.tsv

mapfile -t TASK_IDS < <(
    awk -F '\t' -v chunk="${CHUNK_ID}" 'NR > 1 && $1 == chunk {print $2}' \
        "${CHUNK_MANIFEST}"
)
if [[ ${#TASK_IDS[@]} -eq 0 ]]; then
    echo "No task IDs for chunk ${CHUNK_ID}" >&2
    exit 2
fi

for task_id in "${TASK_IDS[@]}"; do
    "${ENV_PATH}/bin/Rscript" \
        "${H_DIR}/01_estimate_observed_joint_features.R" \
        --run-dir="${RUN_DIR}" --task-id="${task_id}"
done
