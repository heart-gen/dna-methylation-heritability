#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=lgv_geometry_combine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail
RUN_DIR=${LGV_GEOMETRY_DIR:?LGV_GEOMETRY_DIR is required}
H_DIR=${LGV_H_DIR:?LGV_H_DIR is required}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
"${ENV_PATH}/bin/Rscript" "${H_DIR}/13_combine_locus_geometry.R" \
    --run-dir="${RUN_DIR}"
