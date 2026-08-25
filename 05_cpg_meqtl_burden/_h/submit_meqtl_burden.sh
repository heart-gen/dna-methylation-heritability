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
OPEN_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    EXTRA_ARGS="--allow-unlocked"
    OPEN_ARGS="--allow-unlocked"
    if [ -n "${CMB_RUN_ID_OVERRIDE:-}" ]; then
        OPEN_ARGS="$OPEN_ARGS --run-id ${CMB_RUN_ID_OVERRIDE}"
    fi
    log_message "SMOKE RUN: unaccepted upstreams and unlocked keys permitted"
fi

# SMOKE_CHROMS restricts the mapping array to a subset of autosomes (e.g. "22").
# The array bound is otherwise the config autosome policy, so no chromosome can
# be silently skipped.
ARRAY_SPEC="${SMOKE_CHROMS:-1-22}"

RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
    --cohort "$COHORT" --region "$REGION" ${OPEN_ARGS:-$EXTRA_ARGS} | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 00_new_run.R produced no run ID (the upstream gate most likely refused)" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
log_message "run ${RUN_ID}"

mkdir -p "${RUN_DIR}/code"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"

# DRY_RUN inspects the job graph. It stops BEFORE preparing the CpG set,
# because that stage reads every autosome's methylation matrix and takes many
# minutes -- work that tells you nothing about whether the graph is right.
if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1, not submitting. Run dir: ${RUN_DIR}"
    cat <<GRAPH
planned job graph for ${RUN_ID}:
  1   01_prepare_cpg_set.R        submit host, before the array
  2   step_2_map_meqtl.sh         array ${ARRAY_SPEC}, one task per autosome
  2b  02b_combine_meqtl.R         afterany:2    (reconciles cancelled tasks)
  4   04_qc_plots.py              afterok:2b    (lambda, from the nominal pass)
  3   03_vmr_burden.R             afterok:4     (gates on that lambda)
  5   04_check_burden.R           afterok:3
  6   05_finalize_run.R           afterok:5
GRAPH
    exit 0
fi

# The CpG set is prepared on the submit host, not in the array: every mapping
# task reads it, and 22 tasks racing to build it would be both wasteful and
# non-deterministic.
run_r "${SCRIPT_DIR}/01_prepare_cpg_set.R" --run-id "$RUN_ID"

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
EXPORT="ALL,CMB_RUN_ID=${RUN_ID}"

# Stage 2: one array task per autosome, each preparing and mapping its own
# chromosome.
MAP_JOB=$(sbatch --parsable \
    --array="${ARRAY_SPEC}" \
    --chdir="${RUN_DIR}/logs" \
    --export="$EXPORT" \
    "${SCRIPT_DIR}/step_2_map_meqtl.sh")
printf '2\tstep_2_map_meqtl.sh\t%s\n' "$MAP_JOB" >> "$JOBS_TSV"
log_message "mapping array ${MAP_JOB} (${ARRAY_SPEC})"

# Stage 2b: pool the chromosomes and apply the per-region FDR ONCE. afterany,
# so a cancelled chromosome is reconciled rather than leaving the audit blocked.
sbatch_step () {  # name deps cpus mem time command
    sbatch --parsable --dependency="$2" --chdir="${RUN_DIR}/logs" \
        --account=p32505 --partition=short --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%j.out --error=%x-%j.err --export="$EXPORT" \
        --wrap="$6"
}
R_ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"
PY_TQTL="conda run --no-capture-output -p /projects/p32505/opt/envs/genomics python"

COMB_JOB=$(sbatch_step cmb-combine "afterany:${MAP_JOB}" 2 32G 02:00:00 \
    "${R_ENV_SRC} && run_r ${SCRIPT_DIR}/02b_combine_meqtl.R --run-id ${RUN_ID}")
printf '2b\t02b_combine_meqtl.R\t%s\n' "$COMB_JOB" >> "$JOBS_TSV"

# QC runs BEFORE burden, not beside it. AGENTS.md 7.5 makes genomic inflation a
# gate on the burden itself, and lambda can only be estimated from the NOMINAL
# pass -- the permutation table holds one selected top variant per CpG, whose
# p-values are extreme by construction and say nothing about inflation. So
# 04_qc_plots.py computes lambda from the nominal pairs and 03_vmr_burden.R
# reads it back and refuses to aggregate an inflated scan.
QC_JOB=$(sbatch_step cmb-qc "afterok:${COMB_JOB}" 2 32G 02:00:00 \
    "${R_ENV_SRC} && ${PY_TQTL} ${SCRIPT_DIR}/04_qc_plots.py --run-id ${RUN_ID}")
printf '4\t04_qc_plots.py\t%s\n' "$QC_JOB" >> "$JOBS_TSV"

BURDEN_JOB=$(sbatch_step cmb-burden "afterok:${QC_JOB}" 4 32G 02:00:00 \
    "${R_ENV_SRC} && run_r ${SCRIPT_DIR}/03_vmr_burden.R --run-id ${RUN_ID}")
printf '3\t03_vmr_burden.R\t%s\n' "$BURDEN_JOB" >> "$JOBS_TSV"

# The gate reads both the burden table and the genomic-inflation TSV, so it
# waits for both branches.
CHECK_JOB=$(sbatch_step cmb-check "afterok:${BURDEN_JOB}:${QC_JOB}" 1 16G 01:00:00 \
    "${R_ENV_SRC} && run_r ${SCRIPT_DIR}/04_check_burden.R --run-id ${RUN_ID}")
printf '5\t04_check_burden.R\t%s\n' "$CHECK_JOB" >> "$JOBS_TSV"

FINAL_JOB=$(sbatch_step cmb-final "afterok:${CHECK_JOB}" 1 8G 01:00:00 \
    "${R_ENV_SRC} && run_r ${SCRIPT_DIR}/05_finalize_run.R --run-id ${RUN_ID}")
printf '6\t05_finalize_run.R\t%s\n' "$FINAL_JOB" >> "$JOBS_TSV"

log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
