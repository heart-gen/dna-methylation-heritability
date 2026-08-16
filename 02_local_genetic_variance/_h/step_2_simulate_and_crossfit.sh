#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_sim
#SBATCH --output=calibrated-simulation-analysis/_m/logs/%x.%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
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
if [[ ! -f "${SCRIPT_DIR}/01_simulate_and_crossfit.R" ]]; then
    echo "Analysis script is missing from SCRIPT_DIR: ${SCRIPT_DIR}" >&2
    exit 1
fi
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
MANIFEST=${SCENARIO_MANIFEST:-${ANALYSIS_DIR}/_m/config/scenarios.tsv}
OUTPUT_ROOT=${SIMULATION_OUTPUT_ROOT:-${ANALYSIS_DIR}/_m/raw}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/01_simulate_and_crossfit.R" \
    --manifest="${MANIFEST}" \
    --task-id="${SLURM_ARRAY_TASK_ID}" \
    --output-root="${OUTPUT_ROOT}"
