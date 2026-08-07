#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
ENVIRONMENT_FILE=${ANALYSIS_DIR}/environment.yml

if [[ -e "${ENV_PATH}" ]]; then
    echo "Environment path already exists; refusing to modify it: ${ENV_PATH}" >&2
    echo "Set CAL_H2_ENV to a new path or remove the existing environment deliberately." >&2
    exit 1
fi

conda env create --prefix "${ENV_PATH}" --file "${ENVIRONMENT_FILE}"
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/verify_environment.R"

echo "Created and verified calibration environment: ${ENV_PATH}"
