#!/bin/bash
# Donor-robust inference for the DLPFC-hippocampus meQTL-burden slope.
#
#   ./submit_slope_inference.sh AA
#
# Environment:
#   SMOKE_N=1     smoke run (config smoke_bootstrap_n draws; never sealed)
#   DRY_RUN=1     open the run and print the job graph, submit nothing
#
# Chains: open run + write paired draws -> mapping array (chromosome x draw
# block) -> per-draw burden counts -> slope inference, which seals the run.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/05_cpg_meqtl_burden"
COHORT=${1:?usage: submit_slope_inference.sh <cohort>}

OPEN_ARGS=""
if [ -n "${SMOKE_N:-}" ]; then
    OPEN_ARGS="--allow-unlocked"
    [ -n "${CMB_RUN_ID_OVERRIDE:-}" ] && OPEN_ARGS="$OPEN_ARGS --run-id ${CMB_RUN_ID_OVERRIDE}"
    log_message "SMOKE RUN: smoke_bootstrap_n draws, not sealed"
fi

RUN_ID=$(run_r "${SCRIPT_DIR}/06_slope_inference_new_run.R" \
    --cohort "$COHORT" ${OPEN_ARGS} | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 06_slope_inference_new_run.R produced no run ID" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
mkdir -p "${RUN_DIR}/code" "${RUN_DIR}/logs"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"
RUN_CODE="${RUN_DIR}/code/_h"

B=$(awk -F'\t' '$1=="bootstrap_n"{print $2}' "${RUN_DIR}/manifest.tsv")
PER_TASK=$(conda run -p "${V2_ENV_PY}" python -c \
    "import yaml;print(yaml.safe_load(open('${REPO_DIR}/config/meqtl_parameters.yml'))['cross_region_slope_inference']['draws_per_task'])")
N_BLOCKS=$(( (B + PER_TASK - 1) / PER_TASK ))
ARRAY_SPEC="0-$(( 22 * N_BLOCKS - 1 ))"
log_message "run ${RUN_ID}: ${B} draws, ${N_BLOCKS} block(s) of ${PER_TASK}, array ${ARRAY_SPEC}"

if [ "${DRY_RUN:-0}" = "1" ]; then
    cat <<GRAPH
planned job graph for ${RUN_ID}:
  7   step_7_bootstrap_map.sh   array ${ARRAY_SPEC} (chromosome x draw block)
  8a  08a_bootstrap_counts.py   afterok:7
  8   08_slope_inference.R      afterok:8a (seals unless smoke)
GRAPH
    exit 0
fi

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
EXPORT="ALL,CMB_RUN_ID=${RUN_ID},CMB_BOOTSTRAP_N=${B},CMB_DRAWS_PER_TASK=${PER_TASK},V2_RUN_CODE=${RUN_CODE}"

MAP_JOB=$(sbatch --parsable --array="${ARRAY_SPEC}%120" --chdir="${RUN_DIR}/logs" \
    --export="$EXPORT" "${RUN_CODE}/step_7_bootstrap_map.sh")
printf '7\tstep_7_bootstrap_map.sh\t%s\n' "$MAP_JOB" >> "$JOBS_TSV"

sbatch_step () {  # name deps cpus mem time command
    sbatch --parsable --dependency="$2" --chdir="${RUN_DIR}/logs" \
        --account=b1042 --partition=genomics --qos=buyin --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%j.out --error=%x-%j.err --export="$EXPORT" \
        --wrap="$6"
}
R_ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"
PY="conda run --no-capture-output -p ${V2_ENV_PY} python"

COUNT_JOB=$(sbatch_step cmb-bootcount "afterok:${MAP_JOB}" 2 32G 04:00:00 \
    "${R_ENV_SRC} && ${PY} ${RUN_CODE}/08a_bootstrap_counts.py --run-id ${RUN_ID}")
printf '8a\t08a_bootstrap_counts.py\t%s\n' "$COUNT_JOB" >> "$JOBS_TSV"

INF_JOB=$(sbatch_step cmb-slope "afterok:${COUNT_JOB}" 2 16G 02:00:00 \
    "${R_ENV_SRC} && run_r ${RUN_CODE}/08_slope_inference.R --run-id ${RUN_ID}")
printf '8\t08_slope_inference.R\t%s\n' "$INF_JOB" >> "$JOBS_TSV"

log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
