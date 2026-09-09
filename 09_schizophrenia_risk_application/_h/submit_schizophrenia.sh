#!/bin/bash
# Submit one cohort x region cell of 09_schizophrenia_risk_application.
#
#   ./submit_schizophrenia.sh <AA|all_individuals> <region>
#
# Environment:
#   SMOKE_N=1   smoke run (unlocked keys, unaccepted upstreams permitted)
#   DRY_RUN=1   build the job graph, submit nothing
#
# Chain: open run -> define loci -> link VMRs -> [22 risk-variant meQTL tasks]
#        -> pool and correct once -> architecture axis, integration, coupling,
#        GTEx support -> [22 coloc tasks] -> pool coloc -> prioritize ->
#        SuSiE sensitivity -> gate -> panels -> seal.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/09_schizophrenia_risk_application"

COHORT=${1:?usage: submit_schizophrenia.sh <AA|all_individuals> <region>}
REGION=${2:?usage: submit_schizophrenia.sh <AA|all_individuals> <region>}

OPEN_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    OPEN_ARGS="--allow-unlocked"
    [ -n "${SCZ_RUN_ID_OVERRIDE:-}" ] && OPEN_ARGS="$OPEN_ARGS --run-id ${SCZ_RUN_ID_OVERRIDE}"
    log_message "SMOKE RUN: unaccepted upstreams and unlocked keys permitted"
fi

RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
    --cohort "$COHORT" --region "$REGION" ${OPEN_ARGS} | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 00_new_run.R produced no run ID (an upstream gate most likely refused)" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
log_message "run ${RUN_ID}"

mkdir -p "${RUN_DIR}/code"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"
# Snapshot config alongside the code. Stages used to re-read config/ from the
# live working tree, so a config edit -- or a branch switch that removed the
# file -- changed or broke a run already in flight while its manifest attested
# to the original checksum.
mkdir -p "${RUN_DIR}/code/config"
cp -a "${REPO_DIR}/config/." "${RUN_DIR}/code/config/"
RUN_CODE="${RUN_DIR}/code/_h"

if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1, not submitting. Run dir: ${RUN_DIR}"
    cat <<GRAPH
planned job graph for ${RUN_ID}:
  1  01_define_scz_loci.R            submit host (PGC3 loci, index SNPs, hg38 GWAS slice)
  2  02_link_loci_to_vmrs.R          submit host (tested universe on accepted VMRs)
  3  step_3_risk_variant_meqtl.sh    array 1-22: risk variant x CpG meQTL slices
  4  03b_combine_risk_variant_tests.R afterok:3  (one FDR family, applied once)
  5  04_architecture_axis.R          afterok:4
  6  05_integration_enrichment.R     afterok:4
  7  06_transcriptional_coupling.R   afterok:4
  8  07_gtex_support.py              afterok:4
  9  step_8_coloc.sh                 array 1-22, afterok:4 (prepare regions + coloc.abf)
 10  09b_combine_coloc.R             afterok:9
 11  10_prioritize_loci.R            afterok:5,6,7,8,10
 12  11_coloc_susie_sensitivity.R    afterok:11
 13  12_apply_gates.R                afterok:12
 14  13_plot_locus_panels.py         afterok:13
 15  14_finalize_run.R               afterok:14
GRAPH
    exit 0
fi

# Locus definition and VMR linking run on the submit host: minutes of work, and
# every array task reads their output.
log_message "defining PGC3 loci and lifting the GWAS slice to hg38"
run_r "${RUN_CODE}/01_define_scz_loci.R" --run-id "$RUN_ID"
log_message "linking loci to corrected VMRs"
run_r "${RUN_CODE}/02_link_loci_to_vmrs.R" --run-id "$RUN_ID"

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
BASE_EXPORT="ALL,SCZ_RUN_ID=${RUN_ID},V2_REPO_ROOT=${REPO_DIR},V2_RUN_CODE=${RUN_CODE}"

sbatch_step () {  # name deps cpus mem time command
    sbatch --parsable --dependency="$2" --chdir="${RUN_DIR}/logs" \
        --account=b1042 --partition=genomics --qos=buyin --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%j.out --error=%x-%j.err --export="$BASE_EXPORT" \
        --wrap="$6"
}
ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"

RV_JOB=$(sbatch --parsable --chdir="${RUN_DIR}/logs" --export="$BASE_EXPORT" \
    "${RUN_CODE}/step_3_risk_variant_meqtl.sh")
printf '3\tstep_3_risk_variant_meqtl.sh\t%s\n' "$RV_JOB" >> "$JOBS_TSV"

# The FDR family spans all autosomes, so it is corrected once, here, after the
# array completes -- never inside a per-chromosome task.
POOL_JOB=$(sbatch_step scz-pool "afterok:${RV_JOB}" 2 64G 02:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/03b_combine_risk_variant_tests.R --run-id ${RUN_ID}")
printf '4\t03b_combine_risk_variant_tests.R\t%s\n' "$POOL_JOB" >> "$JOBS_TSV"

AXIS_JOB=$(sbatch_step scz-axis "afterok:${POOL_JOB}" 2 32G 02:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/04_architecture_axis.R --run-id ${RUN_ID}")
printf '5\t04_architecture_axis.R\t%s\n' "$AXIS_JOB" >> "$JOBS_TSV"

INTEG_JOB=$(sbatch_step scz-integration "afterok:${POOL_JOB}" 2 32G 02:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/05_integration_enrichment.R --run-id ${RUN_ID}")
printf '6\t05_integration_enrichment.R\t%s\n' "$INTEG_JOB" >> "$JOBS_TSV"

COUP_JOB=$(sbatch_step scz-coupling "afterok:${POOL_JOB}" 2 32G 02:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/06_transcriptional_coupling.R --run-id ${RUN_ID}")
printf '7\t06_transcriptional_coupling.R\t%s\n' "$COUP_JOB" >> "$JOBS_TSV"

GTEX_JOB=$(sbatch_step scz-gtex "afterok:${POOL_JOB}" 4 96G 06:00:00 \
    "${ENV_SRC} && run_py ${RUN_CODE}/07_gtex_support.py --run-id ${RUN_ID}")
printf '8\t07_gtex_support.py\t%s\n' "$GTEX_JOB" >> "$JOBS_TSV"

COLOC_JOB=$(sbatch --parsable --dependency="afterok:${POOL_JOB}" \
    --chdir="${RUN_DIR}/logs" --export="$BASE_EXPORT" \
    "${RUN_CODE}/step_8_coloc.sh")
printf '9\tstep_8_coloc.sh\t%s\n' "$COLOC_JOB" >> "$JOBS_TSV"

COLOC_POOL_JOB=$(sbatch_step scz-coloc-pool "afterok:${COLOC_JOB}" 2 48G 02:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/09b_combine_coloc.R --run-id ${RUN_ID}")
printf '10\t09b_combine_coloc.R\t%s\n' "$COLOC_POOL_JOB" >> "$JOBS_TSV"

PRIO_DEPS="afterok:${AXIS_JOB}:${INTEG_JOB}:${COUP_JOB}:${GTEX_JOB}:${COLOC_POOL_JOB}"
PRIO_JOB=$(sbatch_step scz-prioritize "$PRIO_DEPS" 2 32G 01:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/10_prioritize_loci.R --run-id ${RUN_ID}")
printf '11\t10_prioritize_loci.R\t%s\n' "$PRIO_JOB" >> "$JOBS_TSV"

SUSIE_JOB=$(sbatch_step scz-susie "afterok:${PRIO_JOB}" 2 48G 04:00:00 \
    "${ENV_SRC} && run_r_coloc ${RUN_CODE}/11_coloc_susie_sensitivity.R --run-id ${RUN_ID}")
printf '12\t11_coloc_susie_sensitivity.R\t%s\n' "$SUSIE_JOB" >> "$JOBS_TSV"

GATE_JOB=$(sbatch_step scz-gate "afterok:${SUSIE_JOB}" 2 32G 01:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/12_apply_gates.R --run-id ${RUN_ID}")
printf '13\t12_apply_gates.R\t%s\n' "$GATE_JOB" >> "$JOBS_TSV"

PLOT_JOB=$(sbatch_step scz-panels "afterok:${GATE_JOB}" 2 32G 02:00:00 \
    "${ENV_SRC} && run_py ${RUN_CODE}/13_plot_locus_panels.py --run-id ${RUN_ID}")
printf '14\t13_plot_locus_panels.py\t%s\n' "$PLOT_JOB" >> "$JOBS_TSV"

FINAL_JOB=$(sbatch_step scz-final "afterok:${PLOT_JOB}" 1 8G 01:00:00 \
    "${ENV_SRC} && run_r ${RUN_CODE}/14_finalize_run.R --run-id ${RUN_ID}")
printf '15\t14_finalize_run.R\t%s\n' "$FINAL_JOB" >> "$JOBS_TSV"

log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
