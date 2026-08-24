#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=lgv_regime_manifest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
GEOMETRY_DIR=${LGV_GEOMETRY_DIR:?LGV_GEOMETRY_DIR is required}
H_DIR=${LGV_H_DIR:?LGV_H_DIR is required}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
"${ENV_PATH}/bin/Rscript" "${H_DIR}/07_make_observed_regime_manifest.R" \
    --run-id="${LGV_REGIME_RUN_ID:?}" \
    --features="${GEOMETRY_DIR}/results/combined/locus-geometry.tsv" \
    --cohort="${LGV_COHORT:?}" --region="${LGV_REGION:?}" \
    --vmr-run-id="${LGV_VMR_RUN_ID:?}"
