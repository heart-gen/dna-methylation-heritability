#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
CONFIG=${1:-${ANALYSIS_DIR}/config/analysis.tsv}
OUTPUT=${2:-${ANALYSIS_DIR}/_m/config/scenarios.tsv}

mkdir -p "$(dirname "${OUTPUT}")" "${ANALYSIS_DIR}/_m/logs"
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/00_make_manifest.R" \
    --config="${CONFIG}" \
    --output="${OUTPUT}"
