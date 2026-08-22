#!/bin/bash

## Submit the only active Module 02 production workflow.
## Usage: submit_observed_local_control.sh RUN_ID COHORT REGION VMR_RUN_ID
##
## Stage 00 runs on the submit host, because the chunk count it writes is what
## sizes the Stage 01 array. Every later stage is its own step_*.sh script with
## its own #SBATCH resource block, so a stage can be resubmitted by hand
## without editing this launcher:
##
##   LGV_RUN_DIR=<run dir> LGV_H_DIR=<_h dir> sbatch _h/step_04_derive_score.sh
##
## Environment:
##   LGV_SMOKE_N          restrict the task universe (0 = full run)
##   LGV_VMRS_PER_CHUNK   VMRs per Stage 01 array task (default 5)
##   LGV_MAX_CONCURRENT   Stage 01 array throttle (default 50)
##   DRY_RUN=TRUE         prepare the run, submit nothing

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 RUN_ID COHORT REGION VMR_RUN_ID" >&2
    exit 2
fi
RUN_ID=$1
COHORT=$2
REGION=$3
VMR_RUN_ID=$4

H_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
R_BIN=${ENV_PATH}/bin/Rscript
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
PARTITION=${LGV_PARTITION:-short}
MAX_CONCURRENT=${LGV_MAX_CONCURRENT:-50}
VMRS_PER_CHUNK=${LGV_VMRS_PER_CHUNK:-5}
SMOKE_N=${LGV_SMOKE_N:-0}
DRY_RUN=${DRY_RUN:-FALSE}

if [[ ! -x "${R_BIN}" ]]; then
    echo "Rscript is unavailable: ${R_BIN}" >&2
    exit 1
fi

RUN_DIR=$("${R_BIN}" "${H_DIR}/00_prepare_observed_run.R" \
    --run-id="${RUN_ID}" --cohort="${COHORT}" --region="${REGION}" \
    --vmr-run-id="${VMR_RUN_ID}" --smoke-n="${SMOKE_N}" \
    --vmrs-per-chunk="${VMRS_PER_CHUNK}")
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

submit_step() {
    local dependency=$1
    local step=$2
    shift 2
    local job
    job=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
        ${dependency:+--dependency="${dependency}"} \
        --job-name="lgv_${COHORT}_${REGION}_${step%.sh}" \
        --output="${RUN_DIR}/logs/%x.%j.log" \
        --export="${EXPORTS}" "$@" "${H_DIR}/${step}")
    printf '%s' "${job%%;*}"
}

FEATURE_JOB=$(submit_step "" "step_01_observed_joint_features.sh" \
    --array="1-${N_CHUNKS}%${MAX_CONCURRENT}" \
    --output="${RUN_DIR}/logs/%x.%A_%a.log")

## afterany: a cancelled or failed Stage 01 array must still be reconciled.
COMBINE_JOB=$(submit_step "afterany:${FEATURE_JOB}" "step_02_combine_features.sh")
MODEL_JOB=$(submit_step "afterok:${COMBINE_JOB}" "step_03_apply_joint_model.sh")
SCORE_JOB=$(submit_step "afterok:${MODEL_JOB}" "step_04_derive_score.sh")
QC_JOB=$(submit_step "afterok:${SCORE_JOB}" "step_05_check_score.sh")
FINAL_JOB=$(submit_step "afterok:${QC_JOB}" "step_06_finalize_run.sh")

printf 'stage\tstep_script\tjob_id\n' > "${RUN_DIR}/submitted-jobs.tsv"
printf '%s\t%s\t%s\n' \
    features step_01_observed_joint_features.sh "${FEATURE_JOB}" \
    combine  step_02_combine_features.sh        "${COMBINE_JOB}" \
    model    step_03_apply_joint_model.sh       "${MODEL_JOB}" \
    score    step_04_derive_score.sh            "${SCORE_JOB}" \
    qc       step_05_check_score.sh             "${QC_JOB}" \
    finalize step_06_finalize_run.sh            "${FINAL_JOB}" \
    >> "${RUN_DIR}/submitted-jobs.tsv"

cat <<MSG

Submitted ${RUN_ID} (${COHORT}/${REGION}, ${N_CHUNKS} chunks)

  run dir : ${RUN_DIR}
  jobs    : features=${FEATURE_JOB} combine=${COMBINE_JOB} model=${MODEL_JOB}
            score=${SCORE_JOB} qc=${QC_JOB} finalize=${FINAL_JOB}

A completed job is not acceptance. Stage 05 must return
PASS_RELATIVE_SCORE_OBSERVED_QC and the run must be entered by hand in the
module README's accepted-runs table before any downstream module reads it.
MSG
