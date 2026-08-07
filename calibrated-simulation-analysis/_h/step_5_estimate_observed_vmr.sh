#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_vmr
#SBATCH --output=calibrated-simulation-analysis/_m/logs/%x.%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=04:00:00

set -euo pipefail

if [[ -n "${CAL_H2_SCRIPT_DIR:-}" ]]; then
    SCRIPT_DIR=$(readlink -f "${CAL_H2_SCRIPT_DIR}")
elif [[ -n "${SLURM_SUBMIT_DIR:-}" &&
        -d "${SLURM_SUBMIT_DIR}/calibrated-simulation-analysis/_h" ]]; then
    SCRIPT_DIR=$(readlink -f \
        "${SLURM_SUBMIT_DIR}/calibrated-simulation-analysis/_h")
else
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fi
ANALYSIS_DIR=${CAL_H2_ANALYSIS_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}
if [[ ! -f "${SCRIPT_DIR}/04_estimate_observed_vmr.R" ]]; then
    echo "Analysis script is missing from SCRIPT_DIR: ${SCRIPT_DIR}" >&2
    exit 1
fi
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
REPO_ROOT=${CAL_H2_REPO_ROOT:-$(cd "${ANALYSIS_DIR}/.." && pwd)}
CALIBRATION_MODEL=${CAL_H2_CALIBRATION_MODEL:?CAL_H2_CALIBRATION_MODEL must be set}
OUTPUT_ROOT=${CAL_H2_OBSERVED_OUTPUT_ROOT:?CAL_H2_OBSERVED_OUTPUT_ROOT must be set}
PLINK_ROOT=${CAL_H2_PLINK_ROOT:-/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/all_individuals}
WRITE_DIAGNOSTICS=${CAL_H2_WRITE_DIAGNOSTICS:-FALSE}
: "${REGION:?REGION must be set}"
: "${POPULATION:?POPULATION must be set}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

set +e
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/04_estimate_observed_vmr.R" \
    --region="${REGION}" \
    --population="${POPULATION}" \
    --task-id="${SLURM_ARRAY_TASK_ID}" \
    --repo-root="${REPO_ROOT}" \
    --plink-root="${PLINK_ROOT}" \
    --calibration-model="${CALIBRATION_MODEL}" \
    --output-root="${OUTPUT_ROOT}" \
    --write-diagnostics="${WRITE_DIAGNOSTICS}"
status=$?
set -e
if (( status != 0 )); then
    failure_dir=${OUTPUT_ROOT}/${REGION,,}/${POPULATION}/failures
    mkdir -p "${failure_dir}"
    printf 'task_id\tregion\tpopulation\texit_status\tlog_file\n%s\t%s\t%s\t%s\t%s\n' \
        "${SLURM_ARRAY_TASK_ID}" "${REGION,,}" "${POPULATION}" "${status}" \
        "${SLURM_JOB_NAME:-cal_h2_vmr}.${SLURM_ARRAY_JOB_ID:-NA}_${SLURM_ARRAY_TASK_ID}.log" \
        > "${failure_dir}/vmr-$(printf '%07d' "${SLURM_ARRAY_TASK_ID}").tsv"
    echo "Recorded failed VMR task ${SLURM_ARRAY_TASK_ID} with exit status ${status}" >&2
fi
exit 0
