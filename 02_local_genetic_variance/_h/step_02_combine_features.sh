#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=lgv_combine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#
# Module 02 step 02: reconcile the Stage 01 task universe and combine the
# per-VMR terminal rows.
#
# This step is submitted with `afterany`, not `afterok`: a cancelled or failed
# Stage 01 array must still be reconciled, because the audit record is the
# reconciliation table. Later steps use `afterok` and therefore stop here if
# reconciliation fails.

set -euo pipefail

RUN_DIR=${LGV_RUN_DIR:?LGV_RUN_DIR is required}
H_DIR=${LGV_H_DIR:?LGV_H_DIR is required}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}

"${ENV_PATH}/bin/Rscript" \
    "${H_DIR}/02_combine_observed_joint_features.R" --run-dir="${RUN_DIR}"
