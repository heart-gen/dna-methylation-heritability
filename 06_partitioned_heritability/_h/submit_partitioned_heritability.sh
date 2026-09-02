#!/bin/bash
# Submit one cohort x region cell of 06_partitioned_heritability.
#
#   ./submit_partitioned_heritability.sh <AA|all_individuals> <region>
#
# Environment:
#   SMOKE_N=1        smoke run (unlocked keys, unaccepted upstream permitted)
#   SMOKE_CHROMS=22  restrict the LD-score array (e.g. "22" or "21-22")
#   DRY_RUN=1        build the job graph, submit nothing
#
# Chain: open run -> annotation BED -> liftover -> [LD-score array 1-22]
#        -> munge sumstats -> [S-LDSC per trait] -> FDR/gates -> plot -> seal.
#
# The annotation and liftover run on the submit host: they are seconds of work,
# every array task reads their output, and 22 tasks racing to build the same BED
# would be both wasteful and non-deterministic.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/06_partitioned_heritability"

COHORT=${1:?usage: submit_partitioned_heritability.sh <AA|all_individuals> <region>}
REGION=${2:?usage: submit_partitioned_heritability.sh <AA|all_individuals> <region>}

PY_ENV="/projects/p32505/opt/envs/genomics"
run_py () { conda run --no-capture-output -p "$PY_ENV" python "$@"; }

OPEN_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    OPEN_ARGS="--allow-unlocked"
    if [ -n "${SLDSC_RUN_ID_OVERRIDE:-}" ]; then
        OPEN_ARGS="$OPEN_ARGS --run-id ${SLDSC_RUN_ID_OVERRIDE}"
    fi
    log_message "SMOKE RUN: unaccepted upstream and unlocked keys permitted"
fi

ARRAY_SPEC="${SMOKE_CHROMS:-1-22}"

# The frozen trait list drives the S-LDSC fan-out, so the job graph cannot
# silently disagree with the FDR family the gate will enforce.
TRAITS=$(run_py -c "
import yaml
c = yaml.safe_load(open('${REPO_DIR}/config/partitioned_heritability.yml'))
print(' '.join(t['name'] for t in c['traits']))")
log_message "trait family: ${TRAITS}"

RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" \
    --cohort "$COHORT" --region "$REGION" ${OPEN_ARGS} | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 00_new_run.R produced no run ID (the upstream gate most likely refused)" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
log_message "run ${RUN_ID}"

mkdir -p "${RUN_DIR}/code"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"

# Everything below executes the SNAPSHOT, not the live _h/. Snapshotting and
# then submitting ${SCRIPT_DIR} would let an edit made while jobs are queued
# change what the run executes, while its manifest attests to the older commit.
RUN_CODE="${RUN_DIR}/code/_h"

if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1, not submitting. Run dir: ${RUN_DIR}"
    cat <<GRAPH
planned job graph for ${RUN_ID}:
  1   01_build_annotation.R       submit host (continuous annotation + guardrails)
  2   02_liftover_annotation.py   submit host (hg38 -> hg19)
  3   05_compute_ldscores.sh      array ${ARRAY_SPEC}, annot + LD scores per chrom
  4   04_munge_sumstats.py        afterok:3   (frozen trait family: ${TRAITS})
  5   06_partition_h2.py          afterok:4   one job per trait
  6   07_fdr_and_gates.R          afterok:5   (FDR over the frozen family)
  7   08_plot.py                  afterok:6
  8   09_finalize_run.R           afterok:7
GRAPH
    exit 0
fi

log_message "building the continuous annotation"
run_r "${RUN_CODE}/01_build_annotation.R" --run-id "$RUN_ID"

log_message "lifting the annotation to hg19"
run_py "${RUN_CODE}/02_liftover_annotation.py" --run-id "$RUN_ID"

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
EXPORT="ALL,SLDSC_RUN_ID=${RUN_ID}"

sbatch_step () {  # name deps cpus mem time command
    sbatch --parsable --dependency="$2" --chdir="${RUN_DIR}/logs" \
        --account=b1042 --partition=genomics --qos=buyin --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%j.out --error=%x-%j.err --export="$EXPORT" \
        --wrap="$6"
}
R_ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"
PY_RUN="conda run --no-capture-output -p ${PY_ENV} python"

LD_JOB=$(sbatch --parsable \
    --array="${ARRAY_SPEC}" \
    --chdir="${RUN_DIR}/logs" \
    --export="$EXPORT" \
    "${RUN_CODE}/05_compute_ldscores.sh")
printf '3\t05_compute_ldscores.sh\t%s\n' "$LD_JOB" >> "$JOBS_TSV"
log_message "LD-score array ${LD_JOB} (${ARRAY_SPEC})"

# Munging is independent of the LD scores, but sequencing it after the array
# keeps the cell to one concurrent heavy job and makes a failed chromosome
# cheap to diagnose before an hour of munging runs.
MUNGE_JOB=$(sbatch_step sldsc-munge "afterok:${LD_JOB}" 2 24G 06:00:00 \
    "${R_ENV_SRC} && ${PY_RUN} ${RUN_CODE}/04_munge_sumstats.py --run-id ${RUN_ID}")
printf '4\t04_munge_sumstats.py\t%s\n' "$MUNGE_JOB" >> "$JOBS_TSV"

H2_DEPS=""
for TRAIT in $TRAITS; do
    J=$(sbatch_step "sldsc-h2-${TRAIT}" "afterok:${MUNGE_JOB}" 2 24G 04:00:00 \
        "${R_ENV_SRC} && ${PY_RUN} ${RUN_CODE}/06_partition_h2.py --run-id ${RUN_ID} --trait ${TRAIT}")
    printf '5\t06_partition_h2.py:%s\t%s\n' "$TRAIT" "$J" >> "$JOBS_TSV"
    H2_DEPS="${H2_DEPS}:${J}"
done

# afterok on every trait: the gate refuses a partial family anyway, so failing
# here is clearer than letting stage 07 report the same thing later.
GATE_JOB=$(sbatch_step sldsc-gate "afterok${H2_DEPS}" 1 16G 01:00:00 \
    "${R_ENV_SRC} && run_r ${RUN_CODE}/07_fdr_and_gates.R --run-id ${RUN_ID}")
printf '6\t07_fdr_and_gates.R\t%s\n' "$GATE_JOB" >> "$JOBS_TSV"

PLOT_JOB=$(sbatch_step sldsc-plot "afterok:${GATE_JOB}" 1 16G 01:00:00 \
    "${R_ENV_SRC} && ${PY_RUN} ${RUN_CODE}/08_plot.py --run-id ${RUN_ID}")
printf '7\t08_plot.py\t%s\n' "$PLOT_JOB" >> "$JOBS_TSV"

FINAL_JOB=$(sbatch_step sldsc-final "afterok:${PLOT_JOB}" 1 8G 01:00:00 \
    "${R_ENV_SRC} && run_r ${RUN_CODE}/09_finalize_run.R --run-id ${RUN_ID}")
printf '8\t09_finalize_run.R\t%s\n' "$FINAL_JOB" >> "$JOBS_TSV"

log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
