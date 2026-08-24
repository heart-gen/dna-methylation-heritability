#!/bin/bash
# Evaluate a completed independent validation run against the frozen model and
# write the terminal decision into a distinct immutable derived run.

#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=joint_pve_decide
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00

set -euo pipefail

ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
SCRIPT_DIR=${JOINT_PVE_SCRIPT_DIR:?}
VALIDATION_ROOT=${JOINT_PVE_VALIDATION_ROOT:?}
DECISION_ROOT=${JOINT_PVE_DECISION_ROOT:?}
TRAIN_ROOT=${JOINT_PVE_TRAIN_ROOT:?}

"${ENV_PATH}/bin/Rscript" \
    "${SCRIPT_DIR}/24_evaluate_joint_pve_validation.R" \
    --input_dir="${VALIDATION_ROOT}/raw" \
    --manifest="${VALIDATION_ROOT}/config/scenarios.tsv" \
    --config="${DECISION_ROOT}/config/joint-pve-effective-recovery.tsv" \
    --criteria="${DECISION_ROOT}/config/joint-pve-acceptance-criteria.tsv" \
    --model="${TRAIN_ROOT}/combined/joint-pve-calibrator.rds" \
    --model_sha256="${TRAIN_ROOT}/combined/joint-pve-calibrator.sha256" \
    --output_dir="${DECISION_ROOT}/combined" \
    --fail_on_rejection=FALSE

find "${DECISION_ROOT}/combined" -maxdepth 1 -type f -print0 | \
    sort -z | xargs -0 sha256sum \
    > "${DECISION_ROOT}/provenance/output-checksums.sha256"
