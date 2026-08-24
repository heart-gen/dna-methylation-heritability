#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=lgv_regime
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

RUN_DIR=${LGV_RUN_DIR:?LGV_RUN_DIR is required}
H_DIR=${LGV_H_DIR:?LGV_H_DIR is required}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
CHUNK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}

"${ENV_PATH}/bin/Rscript" \
    "${H_DIR}/08_run_observed_regime_scenario.R" \
    --run-dir="${RUN_DIR}" --chunk-id="${CHUNK_ID}"
