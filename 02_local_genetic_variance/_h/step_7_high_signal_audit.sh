#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_audit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail

SCRIPT_DIR=${CAL_H2_SCRIPT_DIR:?CAL_H2_SCRIPT_DIR must be set}
REPO_ROOT=${CAL_H2_REPO_ROOT:?CAL_H2_REPO_ROOT must be set}
MANIFEST=${CAL_H2_AUDIT_MANIFEST:?CAL_H2_AUDIT_MANIFEST must be set}
OUTPUT_ROOT=${CAL_H2_AUDIT_OUTPUT_ROOT:?CAL_H2_AUDIT_OUTPUT_ROOT must be set}
CALIBRATION_MODEL=${CAL_H2_CALIBRATION_MODEL:?CAL_H2_CALIBRATION_MODEL must be set}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
AUDIT_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set}

record=$(awk -F '\t' -v row="${AUDIT_ID}" 'NR == row + 1 {print; exit}' "${MANIFEST}")
[[ -n "${record}" ]] || { echo "Audit ID is outside manifest: ${AUDIT_ID}" >&2; exit 1; }
IFS=$'\t' read -r audit_id population region task_id vmr_id upstream_run \
    vmr_set_id original_r2 original_rho2 sensitivity excluded_fids \
    seed_repeat fold_seed <<< "${record}"
[[ "${audit_id}" == "${AUDIT_ID}" ]] || {
    echo "Manifest row/audit ID mismatch" >&2; exit 1; }

VMR_RUN_DIR=${REPO_ROOT}/01_vmr_catalog/_m/runs/${upstream_run}
TASK_OUTPUT=${OUTPUT_ROOT}/tasks/audit-$(printf '%05d' "${AUDIT_ID}")
mkdir -p "${OUTPUT_ROOT}/summaries"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/04_estimate_observed_vmr.R" \
    --region="${region}" \
    --population="${population}" \
    --cohort="${population}" \
    --task-id="${task_id}" \
    --repo-root="${REPO_ROOT}" \
    --vmr-run-dir="${VMR_RUN_DIR}" \
    --calibration-model="${CALIBRATION_MODEL}" \
    --output-root="${TASK_OUTPUT}" \
    --write-diagnostics=TRUE \
    --fold-seed="${fold_seed}" \
    --exclude-fids="${excluded_fids}"

summary=${TASK_OUTPUT}/${region}/${population}/summary/vmr-$(printf '%07d' "${task_id}").tsv
[[ -s "${summary}" ]] || { echo "Audit summary was not created: ${summary}" >&2; exit 1; }
cp "${summary}" "${OUTPUT_ROOT}/summaries/audit-$(printf '%05d' "${AUDIT_ID}").tsv"
