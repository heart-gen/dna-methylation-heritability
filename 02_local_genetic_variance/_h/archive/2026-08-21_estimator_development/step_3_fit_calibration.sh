#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_fit
#SBATCH --output=calibrated-simulation-analysis/_m/logs/%x.%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00

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
if [[ ! -f "${SCRIPT_DIR}/02_fit_calibration.R" ]]; then
    echo "Analysis script is missing from SCRIPT_DIR: ${SCRIPT_DIR}" >&2
    exit 1
fi
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
RUN_ROOT=${CAL_H2_RUN_ROOT:-${ANALYSIS_DIR}/_m}

"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/02_fit_calibration.R" \
    --input="${RUN_ROOT}/raw/calibration" \
    --criteria="${RUN_ROOT}/config/acceptance-criteria.tsv" \
    --analysis-config="${RUN_ROOT}/config/analysis.tsv" \
    --output-model="${RUN_ROOT}/calibration/elastic-net-calibration.rds" \
    --output-manifest="${RUN_ROOT}/calibration/calibration-manifest.tsv" \
    --output-tuning="${RUN_ROOT}/calibration/hybrid-weight-tuning.tsv" \
    --session-info="${RUN_ROOT}/calibration/session-info.txt"
