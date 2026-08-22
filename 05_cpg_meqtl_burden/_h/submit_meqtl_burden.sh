#!/bin/bash
# Submit one cohort x region cell of 05_cpg_meqtl_burden.
#
#   ./submit_meqtl_burden.sh <AA|all_individuals> <caudate|dlpfc|hippocampus>
#
# Environment:
#   SMOKE_N=1     smoke run (unlocked keys, unaccepted upstreams permitted)
#   DRY_RUN=1     build everything, submit nothing
#
# Chains: open run -> prepare CpG set -> per-chromosome mapping array ->
# aggregate to VMR burden. The final two steps run under sbatch dependencies so
# the whole cell is one submission.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/05_cpg_meqtl_burden"

COHORT=${1:?usage: submit_meqtl_burden.sh <AA|all_individuals> <region>}
REGION=${2:?usage: submit_meqtl_burden.sh <AA|all_individuals> <region>}

EXTRA_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    EXTRA_ARGS="--allow-unlocked"
    log_message "SMOKE RUN: unaccepted upstreams and unlocked keys permitted"
fi

RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
    --cohort "$COHORT" --region "$REGION" $EXTRA_ARGS | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 00_new_run.R produced no run ID (the upstream gate most likely refused)" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
log_message "run ${RUN_ID}"

mkdir -p "${RUN_DIR}/code"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"

# The CpG set is prepared on the submit host, not in the array: every mapping
# task reads it, and 22 tasks racing to build it would be both wasteful and
# non-deterministic.
run_r "${SCRIPT_DIR}/01_prepare_cpg_set.R" --run-id "$RUN_ID"

if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1, not submitting. Run dir: ${RUN_DIR}"
    exit 0
fi

MAP_JOB=$(sbatch --parsable \
    --chdir="${RUN_DIR}/logs" \
    --export=ALL,CMB_RUN_ID="$RUN_ID" \
    "${SCRIPT_DIR}/step_2_map_meqtl.sh")
log_message "mapping array ${MAP_JOB}"

BURDEN_JOB=$(sbatch --parsable \
    --dependency=afterok:"$MAP_JOB" \
    --chdir="${RUN_DIR}/logs" \
    --account=b1042 --partition=genomicsguest \
    --job-name=cmb-burden --cpus-per-task=4 --mem=32G --time=02:00:00 \
    --output=%x-%j.out --error=%x-%j.err \
    --wrap="source ${REPO_DIR}/00_shared/slurm.sh && \
            run_r ${SCRIPT_DIR}/03_vmr_burden.R --run-id ${RUN_ID}")
log_message "burden job ${BURDEN_JOB} (afterok:${MAP_JOB})"
