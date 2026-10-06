#!/bin/bash
# Submit one region of 02c_cis_greml_sensitivity.
#
#   ./submit_cis_greml.sh AA <region>
#
# Environment:
#   SMOKE=1                    smoke run (unlocked config, a spread of eligible VMRs)
#   CGS_RUN_ID_OVERRIDE=<id>   smoke only: choose the run ID
#   DRY_RUN=1                  open the run and print the plan only
#
# Chain: open (submit host) -> per-VMR cis-GREML (array) -> summarize -> gate
#        -> plot -> finalize. After all three regions are accepted:
#   Rscript _h/06_collate.R

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/02c_cis_greml_sensitivity"
COHORT=${1:?usage: submit_cis_greml.sh AA <region>}
REGION=${2:?usage: submit_cis_greml.sh AA <region>}

OPEN_ARGS=(--cohort "$COHORT" --region "$REGION")
if [ -n "${SMOKE:-}" ]; then
    OPEN_ARGS+=(--allow-unlocked)
    [ -n "${CGS_RUN_ID_OVERRIDE:-}" ] && OPEN_ARGS+=(--run-id "$CGS_RUN_ID_OVERRIDE")
    log_message "SMOKE RUN: unlocked config permitted"
fi
RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" "${OPEN_ARGS[@]}" | tail -n 1)
[ -z "$RUN_ID" ] && { echo "ERROR: 00_new_run.R produced no run ID" >&2; exit 1; }
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
log_message "run ${RUN_ID}"

mkdir -p "${RUN_DIR}/code/config"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"
cp -a "${REPO_DIR}/config/." "${RUN_DIR}/code/config/"
RUN_CODE="${RUN_DIR}/code/_h"
N_TASKS=$(awk -F'\t' '$1=="n_tasks"{print $2}' "${RUN_DIR}/manifest.tsv")

if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1: ${N_TASKS} fit tasks; not submitting. Run dir: ${RUN_DIR}"
    exit 0
fi

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
BASE_EXPORT="ALL,V2_REPO_ROOT=${REPO_DIR},V2_RUN_CODE=${RUN_CODE}"
ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"
sbatch_step () {  # name deps cpus mem time array command
    local extra=()
    [ -n "$2" ] && extra+=(--dependency="$2")
    [ -n "$6" ] && extra+=(--array="$6")
    sbatch --parsable "${extra[@]}" --chdir="${RUN_DIR}/logs" \
        --account=b1042 --partition=genomics --qos=buyin --job-name="$1" \
        --cpus-per-task="$3" --mem="$4" --time="$5" \
        --output=%x-%A_%a.out --error=%x-%A_%a.err --export="$BASE_EXPORT" \
        --wrap="$7"
}
r_cmd () { echo "${ENV_SRC} && run_r ${RUN_CODE}/$1 --run-id ${RUN_ID}"; }

FIT=$(sbatch_step cgs-fit "" 1 8G 04:00:00 "1-${N_TASKS}" "$(r_cmd 01_fit_cis_greml.R)")
printf '1\t01_fit_cis_greml.R\t%s\n' "$FIT" >> "$JOBS_TSV"
SUM=$(sbatch_step cgs-summary "afterok:${FIT}" 1 16G 01:00:00 "" "$(r_cmd 02_summarize.R)")
printf '2\t02_summarize.R\t%s\n' "$SUM" >> "$JOBS_TSV"
GATE=$(sbatch_step cgs-gate "afterok:${SUM}" 1 4G 00:30:00 "" "$(r_cmd 03_apply_gate.R)")
printf '3\t03_apply_gate.R\t%s\n' "$GATE" >> "$JOBS_TSV"
PLOT=$(sbatch_step cgs-plot "afterok:${GATE}" 1 8G 00:30:00 "" "$(r_cmd 04_plot.R)")
printf '4\t04_plot.R\t%s\n' "$PLOT" >> "$JOBS_TSV"
FIN=$(sbatch_step cgs-final "afterok:${PLOT}" 1 4G 00:30:00 "" "$(r_cmd 05_finalize_run.R)")
printf '5\t05_finalize_run.R\t%s\n' "$FIN" >> "$JOBS_TSV"
log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
