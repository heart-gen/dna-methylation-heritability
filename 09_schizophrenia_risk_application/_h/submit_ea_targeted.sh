#!/bin/bash
# Submit the targeted EA check of the illustrative loci (stages 19-22).
#
#   ./submit_ea_targeted.sh
#
# Environment:
#   SMOKE_N=1   smoke run (unlocked keys permitted; never sealed)
#   DRY_RUN=1   print the job graph, open nothing
#
# Chain: open run + targets (submit host) -> [stage 20 per region]
#        -> stage 21 coloc (coloc env) -> stage 22 read + seal.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/09_schizophrenia_risk_application"
PY_ENV="${V2_ENV_PY}"

OPEN_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    OPEN_ARGS="--allow-unlocked"
    [ -n "${SCZ_RUN_ID_OVERRIDE:-}" ] && OPEN_ARGS="$OPEN_ARGS --run-id ${SCZ_RUN_ID_OVERRIDE}"
    log_message "SMOKE RUN: not sealed"
fi
if [ "${DRY_RUN:-0}" = "1" ]; then
    cat <<GRAPH
planned job graph (targeted EA check):
  19  19_ea_targeted_new_run.R   submit host (targets from the accepted AA runs)
  20  20_ea_targeted_test.py     one job per region
  21  21_ea_coloc.R              afterok:20 (coloc env, non-gating)
  22  22_ea_targeted_summarize.R afterok:21, read + seal
GRAPH
    exit 0
fi

RUN_ID=$(run_r "${SCRIPT_DIR}/19_ea_targeted_new_run.R" ${OPEN_ARGS} | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 19_ea_targeted_new_run.R produced no run ID" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"
RUN_CODE="${RUN_DIR}/code/_h"
log_message "run ${RUN_ID}"

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
EXPORT="ALL,V2_REPO_ROOT=${REPO_DIR},V2_RUN_CODE=${RUN_CODE}"
sbatch_step () {  # name deps cpus mem time command
    sbatch --parsable ${2:+--dependency="$2"} --chdir="${RUN_DIR}/logs" \
        --account=b1042 --partition=genomics --qos=buyin --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%j.out --error=%x-%j.err --export="$EXPORT" \
        --wrap="$6"
}
SRC="source ${REPO_DIR}/00_shared/slurm.sh"

DEPS=""
for REGION in caudate dlpfc hippocampus; do
    J=$(sbatch_step "scz-ea-${REGION}" "" 4 48G 06:00:00 \
        "${SRC} && run_py ${RUN_CODE}/20_ea_targeted_test.py --run-id ${RUN_ID} --region ${REGION} --threads 8")
    printf '20\t20_ea_targeted_test.py:%s\t%s\n' "$REGION" "$J" >> "$JOBS_TSV"
    DEPS="${DEPS}:${J}"
done
COLOC=$(sbatch_step scz-ea-coloc "afterok${DEPS}" 1 16G 01:00:00 \
    "${SRC} && run_r_coloc ${RUN_CODE}/21_ea_coloc.R --run-id ${RUN_ID}")
printf '21\t21_ea_coloc.R\t%s\n' "$COLOC" >> "$JOBS_TSV"
FINAL=$(sbatch_step scz-ea-final "afterok:${COLOC}" 1 8G 01:00:00 \
    "${SRC} && run_r ${RUN_CODE}/22_ea_targeted_summarize.R --run-id ${RUN_ID}")
printf '22\t22_ea_targeted_summarize.R\t%s\n' "$FINAL" >> "$JOBS_TSV"
log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
