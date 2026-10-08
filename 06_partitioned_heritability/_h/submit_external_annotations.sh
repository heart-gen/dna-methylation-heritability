#!/bin/bash
# Submit the external-annotation S-LDSC benchmark (stages 10-13).
#
#   ./submit_external_annotations.sh AA
#
# Environment:
#   SMOKE_N=1        smoke run (unlocked keys permitted; never sealed)
#   SMOKE_CHROMS=22  restrict the LD-score array; stage 12 then fails by design
#                    (it needs all 22), so a chromosome smoke stops after 11
#   SMOKE_TRAITS="scz cad"  restrict the regression fan-out (smoke only)
#   DRY_RUN=1        open nothing, print the job graph
#
# Chain: open run + stage annotations and sumstats (submit host)
#        -> [LD-score array, chromosome x every annotation]
#        -> [S-LDSC, annotation x trait] -> combine + seal.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/06_partitioned_heritability"
COHORT=${1:?usage: submit_external_annotations.sh <cohort>}
PY_ENV="/projects/p32505/opt/envs/genomics"
run_py () { conda run --no-capture-output -p "$PY_ENV" python "$@"; }

OPEN_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    OPEN_ARGS="--allow-unlocked"
    [ -n "${SLDSC_RUN_ID_OVERRIDE:-}" ] && OPEN_ARGS="$OPEN_ARGS --run-id ${SLDSC_RUN_ID_OVERRIDE}"
    log_message "SMOKE RUN: not sealed"
fi
ARRAY_SPEC="${SMOKE_CHROMS:-1-22}"
TRAITS=$(run_py -c "
import yaml
c = yaml.safe_load(open('${REPO_DIR}/config/partitioned_heritability.yml'))
print(' '.join(t['name'] for t in c['traits']))")
if [ -n "${SMOKE_N:-}" ] && [ -n "${SMOKE_TRAITS:-}" ]; then TRAITS="$SMOKE_TRAITS"; fi

if [ "${DRY_RUN:-0}" = "1" ]; then
    cat <<GRAPH
planned job graph (external-annotation benchmark, cohort ${COHORT}):
  10  10_external_new_run.R        submit host (stage annotations + sumstats)
  11  11_external_ldscores.sh      array ${ARRAY_SPEC}, every ready annotation per chromosome
  12  12_external_partition_h2.py  afterok:11, one job per annotation x trait (${TRAITS})
  13  13_external_summarize.R      afterok:12, combine + seal
GRAPH
    exit 0
fi

RUN_ID=$(run_r "${SCRIPT_DIR}/10_external_new_run.R" --cohort "$COHORT" ${OPEN_ARGS} | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 10_external_new_run.R produced no run ID" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"
RUN_CODE="${RUN_DIR}/code/_h"
log_message "run ${RUN_ID}"

ANNOTS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} $h["status"]=="ready"{print $h["annotation"]}' \
    "${RUN_DIR}/annotation/annotation-list.tsv")

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
EXPORT="ALL,SLDSC_RUN_ID=${RUN_ID},V2_REPO_ROOT=${REPO_DIR},V2_RUN_CODE=${RUN_CODE}"
sbatch_step () {  # name deps cpus mem time command
    sbatch --parsable --dependency="$2" --chdir="${RUN_DIR}/logs" \
        --account=b1042 --partition=genomics --qos=buyin --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%j.out --error=%x-%j.err --export="$EXPORT" \
        --wrap="$6"
}
R_ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"
PY_RUN="conda run --no-capture-output -p ${PY_ENV} python"

LD_JOB=$(sbatch --parsable --array="${ARRAY_SPEC}" --chdir="${RUN_DIR}/logs" \
    --export="$EXPORT" "${RUN_CODE}/11_external_ldscores.sh")
printf '11\t11_external_ldscores.sh\t%s\n' "$LD_JOB" >> "$JOBS_TSV"

H2_DEPS=""
for A in $ANNOTS; do
    for T in $TRAITS; do
        J=$(sbatch_step "sldsc-x-${T}" "afterok:${LD_JOB}" 2 24G 04:00:00 \
            "${R_ENV_SRC} && ${PY_RUN} ${RUN_CODE}/12_external_partition_h2.py --run-id ${RUN_ID} --annotation ${A} --trait ${T}")
        printf '12\t12_external_partition_h2.py:%s:%s\t%s\n' "$A" "$T" "$J" >> "$JOBS_TSV"
        H2_DEPS="${H2_DEPS}:${J}"
    done
done

FINAL=$(sbatch_step sldsc-x-final "afterok${H2_DEPS}" 1 8G 01:00:00 \
    "${R_ENV_SRC} && run_r ${RUN_CODE}/13_external_summarize.R --run-id ${RUN_ID}")
printf '13\t13_external_summarize.R\t%s\n' "$FINAL" >> "$JOBS_TSV"
log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
