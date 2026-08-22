#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=lgv_regime_summary
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail

RUN_DIR=${LGV_RUN_DIR:?LGV_RUN_DIR is required}
H_DIR=${LGV_H_DIR:?LGV_H_DIR is required}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}

"${ENV_PATH}/bin/Rscript" \
    "${H_DIR}/09_summarize_observed_regime.R" --run-dir="${RUN_DIR}"
