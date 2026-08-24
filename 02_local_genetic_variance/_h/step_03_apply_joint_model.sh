#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=lgv_model
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#
# Module 02 step 03: verify the frozen joint-model checksum and apply it once.

set -euo pipefail

RUN_DIR=${LGV_RUN_DIR:?LGV_RUN_DIR is required}
H_DIR=${LGV_H_DIR:?LGV_H_DIR is required}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}

"${ENV_PATH}/bin/Rscript" "${H_DIR}/03_apply_frozen_joint_model.R" --run-dir="${RUN_DIR}"
