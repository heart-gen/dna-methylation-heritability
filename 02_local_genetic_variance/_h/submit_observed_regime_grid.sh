#!/bin/bash

## Submit the observed-regime diagnostic grid.
##
## Usage: submit_observed_regime_grid.sh [RUN_ID]
##
## This grid does not change the frozen joint model, its calibration, or its
## acceptance gates. It characterises the existing model on real cis-window
## genotypes at the observed sample size, which the AR(1) validation grid never
## covered, and answers whether the observed lower-boundary mass is an expected
## resolution limit or a simulation-to-data mismatch.
##
## Environment:
##   LGV_MAX_CONCURRENT  array throttle (default 50)
##   DRY_RUN=TRUE        prepare the run, submit nothing

set -euo pipefail

RUN_ID=${1:-}
H_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
R_BIN=${ENV_PATH}/bin/Rscript
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
PARTITION=${LGV_PARTITION:-short}
MAX_CONCURRENT=${LGV_MAX_CONCURRENT:-50}
DRY_RUN=${DRY_RUN:-FALSE}

if [[ ! -x "${R_BIN}" ]]; then
    echo "Rscript is unavailable: ${R_BIN}" >&2
    exit 1
fi

RUN_DIR=$("${R_BIN}" "${H_DIR}/07_make_observed_regime_manifest.R" \
    ${RUN_ID:+--run-id="${RUN_ID}"})
RUN_DIR=$(echo "${RUN_DIR}" | tail -1 | tr -d '[:space:]')
N_CHUNKS=$(awk -F '\t' 'NR > 1 {if ($1 > n) n=$1} END {print n+0}' \
    "${RUN_DIR}/config/chunk-manifest.tsv")
if [[ ${N_CHUNKS} -lt 1 ]]; then
    echo "Prepared run has no chunks" >&2
    exit 1
fi

if [[ "${DRY_RUN,,}" == "true" ]]; then
    echo "Prepared dry run: ${RUN_DIR} (${N_CHUNKS} chunks)"
    exit 0
fi

EXPORTS="ALL,LGV_RUN_DIR=${RUN_DIR},LGV_H_DIR=${H_DIR},CAL_H2_ENV=${ENV_PATH}"

SCENARIO_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --partition="${PARTITION}" \
    --array="1-${N_CHUNKS}%${MAX_CONCURRENT}" \
    --output="${RUN_DIR}/logs/%x.%A_%a.log" \
    --export="${EXPORTS}" \
    "${H_DIR}/step_07_observed_regime_scenarios.sh")
SCENARIO_JOB=${SCENARIO_JOB%%;*}

## afterany: a cancelled or failed array must still be reconciled.
SUMMARY_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --partition="${PARTITION}" \
    --dependency="afterany:${SCENARIO_JOB}" \
    --output="${RUN_DIR}/logs/%x.%j.log" \
    --export="${EXPORTS}" \
    "${H_DIR}/step_08_summarize_observed_regime.sh")
SUMMARY_JOB=${SUMMARY_JOB%%;*}

printf 'stage\tstep_script\tjob_id\n' > "${RUN_DIR}/submitted-jobs.tsv"
printf '%s\t%s\t%s\n' \
    scenarios step_07_observed_regime_scenarios.sh "${SCENARIO_JOB}" \
    summary   step_08_summarize_observed_regime.sh "${SUMMARY_JOB}" \
    >> "${RUN_DIR}/submitted-jobs.tsv"

cat <<MSG

Submitted observed-regime grid (${N_CHUNKS} chunks)

  run dir : ${RUN_DIR}
  jobs    : scenarios=${SCENARIO_JOB} summary=${SUMMARY_JOB}

This grid is diagnostic. It authorises no scoring change on its own; read
results/combined/regime-verdict.tsv before touching Stage 04.
MSG
