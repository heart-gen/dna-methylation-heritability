#!/bin/bash

set -euo pipefail

TEST_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${TEST_DIR}/.." && pwd)
SCRIPT_DIR=${ANALYSIS_DIR}/_h
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
TEST_ROOT=$(mktemp -d /tmp/cal-h2-slurm-spool.XXXXXX)
trap 'rm -rf "${TEST_ROOT}"' EXIT

cp "${SCRIPT_DIR}/step_2_simulate_and_crossfit.sh" \
    "${TEST_ROOT}/slurm_spooled_script"
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/00_make_manifest.R" \
    --config="${ANALYSIS_DIR}/config/smoke.tsv" \
    --output="${TEST_ROOT}/scenarios.tsv"

env \
    CAL_H2_SCRIPT_DIR="${SCRIPT_DIR}" \
    CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}" \
    CAL_H2_ENV="${ENV_PATH}" \
    SCENARIO_MANIFEST="${TEST_ROOT}/scenarios.tsv" \
    SIMULATION_OUTPUT_ROOT="${TEST_ROOT}/raw" \
    SLURM_ARRAY_TASK_ID=1 \
    SLURM_CPUS_PER_TASK=1 \
    bash "${TEST_ROOT}/slurm_spooled_script"

test -s "${TEST_ROOT}/raw/calibration/scenario-0000001.tsv"
echo "SLURM spool-path regression test passed"
