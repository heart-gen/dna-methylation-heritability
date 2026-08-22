#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=lgv_score
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#
# Module 02 step 04: derive the within-cell relative local-genetic-control score.

set -euo pipefail

RUN_DIR=${LGV_RUN_DIR:?LGV_RUN_DIR is required}
H_DIR=${LGV_H_DIR:?LGV_H_DIR is required}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}

"${ENV_PATH}/bin/Rscript" "${H_DIR}/04_derive_local_snp_contribution_score.R" --run-dir="${RUN_DIR}"
