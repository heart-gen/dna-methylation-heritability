#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_fit_boundary
#SBATCH --output=02_local_genetic_variance/_m/logs/%x.%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

if [[ -n "${CAL_H2_SCRIPT_DIR:-}" ]]; then
    SCRIPT_DIR=$(readlink -f "${CAL_H2_SCRIPT_DIR}")
else
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fi
ANALYSIS_DIR=${CAL_H2_ANALYSIS_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
RUN_ROOT=${CAL_H2_RUN_ROOT:-${ANALYSIS_DIR}/_m}
CANDIDATES=${CAL_H2_CANDIDATES:-${RUN_ROOT}/config/boundary-aware-candidates.tsv}

if [[ ! -f "${SCRIPT_DIR}/18_fit_boundary_aware_calibration.R" ]]; then
    echo "Boundary-aware fit script missing from SCRIPT_DIR: ${SCRIPT_DIR}" >&2
    exit 1
fi

"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/18_fit_boundary_aware_calibration.R" \
    --input="${RUN_ROOT}/raw/calibration" \
    --candidates="${CANDIDATES}" \
    --criteria="${RUN_ROOT}/config/acceptance-criteria.tsv" \
    --analysis-config="${RUN_ROOT}/config/analysis.tsv" \
    --output-model="${RUN_ROOT}/calibration/elastic-net-calibration.rds" \
    --output-manifest="${RUN_ROOT}/calibration/calibration-manifest.tsv" \
    --output-tuning="${RUN_ROOT}/calibration/candidate-selection.tsv" \
    --output-selected="${RUN_ROOT}/calibration/selected-candidate.tsv" \
    --session-info="${RUN_ROOT}/calibration/session-info.txt" \
    --fail-closed=TRUE
