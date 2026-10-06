#!/bin/bash
# Submit one run of 02b_greml_simulation_benchmark.
#
#   ./submit_greml_benchmark.sh ar1                 # arm 1, all simulated n
#   ./submit_greml_benchmark.sh observed AA <region> # arm 2, one region
#
# Environment:
#   SMOKE=1                      smoke run (unlocked config, reduced grid)
#   GREML_RUN_ID_OVERRIDE=<id>   smoke only: choose the run ID
#   DRY_RUN=1                    open the run and print the job graph only
#
# Arm 1 chain: inputs + LD scores (array) -> stratified GRMs (array)
#              -> REML (array) -> summarize -> gate -> plot -> finalize
# Arm 2 chain: REML (array over locus chunks) -> summarize -> gate -> plot
#              -> finalize
#
# Arm 1's n = 5,000 and 10,000 cells need far more memory and time than the
# rest, so the GRM and REML arrays are split into a small and a large class
# and submitted with separate resources.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_ROOT="${REPO_DIR}/02b_greml_simulation_benchmark"

ARM=${1:?usage: submit_greml_benchmark.sh ar1 | observed AA <region>}
OPEN_ARGS=(--arm "$ARM")
if [ "$ARM" = "observed" ]; then
    COHORT=${2:?observed needs a cohort}
    REGION=${3:?observed needs a region}
    OPEN_ARGS+=(--cohort "$COHORT" --region "$REGION")
fi
if [ -n "${SMOKE:-}" ]; then
    OPEN_ARGS+=(--allow-unlocked)
    [ -n "${GREML_RUN_ID_OVERRIDE:-}" ] && OPEN_ARGS+=(--run-id "$GREML_RUN_ID_OVERRIDE")
    log_message "SMOKE RUN: unlocked config permitted, reduced grid"
fi

RUN_ID=$(run_r "${SCRIPT_DIR}/00_new_run.R" "${OPEN_ARGS[@]}" | tail -n 1)
if [ -z "$RUN_ID" ]; then
    echo "ERROR: 00_new_run.R produced no run ID" >&2
    exit 1
fi
RUN_DIR="${MODULE_ROOT}/_m/runs/${RUN_ID}"
log_message "run ${RUN_ID}"

# Snapshot code and config so an edit or branch switch cannot change a run in
# flight (see 09b_aging_application/_h/run_config.R).
mkdir -p "${RUN_DIR}/code/config"
cp -a "$SCRIPT_DIR" "${RUN_DIR}/code/_h"
cp -a "${REPO_DIR}/config/." "${RUN_DIR}/code/config/"
RUN_CODE="${RUN_DIR}/code/_h"

# Task-index ranges from a task manifest, as a compact sbatch --array list.
# $1 = manifest, $2 = awk condition on column N (or "1" for all rows).
task_ranges () {
    awk -F'\t' -v cond="$2" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        { N = (("N" in col) ? $col["N"] : 0); t = $col["task"]
          keep = (cond == "1") || (cond == "small" && N <= 1000) || (cond == "large" && N > 1000)
          if (keep) print t }' "$1" | sort -n | awk '
        NR == 1 { s = $1; p = $1; next }
        $1 == p + 1 { p = $1; next }
        { printf "%s%s", sep, (s == p ? s : s "-" p); sep = ","; s = $1; p = $1 }
        END { if (NR) printf "%s%s\n", sep, (s == p ? s : s "-" p) }'
}

JOBS_TSV="${RUN_DIR}/submitted-jobs.tsv"
printf 'step\tscript\tjob_id\n' > "$JOBS_TSV"
BASE_EXPORT="ALL,V2_REPO_ROOT=${REPO_DIR},V2_RUN_CODE=${RUN_CODE}"
ENV_SRC="source ${REPO_DIR}/00_shared/slurm.sh"
record () { printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$JOBS_TSV"; }

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

if [ "${DRY_RUN:-0}" = "1" ]; then
    log_message "DRY_RUN=1, not submitting. Run dir: ${RUN_DIR}"
    if [ "$ARM" = "ar1" ]; then
        echo "ldscore tasks:       $(task_ranges "${RUN_DIR}/config/ldscore-tasks.tsv" 1)"
        echo "stratify small/large: $(task_ranges "${RUN_DIR}/config/stratify-tasks.tsv" small) / $(task_ranges "${RUN_DIR}/config/stratify-tasks.tsv" large)"
        echo "reml small/large:     $(task_ranges "${RUN_DIR}/config/reml-tasks.tsv" small) / $(task_ranges "${RUN_DIR}/config/reml-tasks.tsv" large)"
    else
        echo "reml tasks: $(task_ranges "${RUN_DIR}/config/locus-tasks.tsv" 1)"
    fi
    exit 0
fi

DEPS=""
if [ "$ARM" = "ar1" ]; then
    IN=$(sbatch_step greml-inputs "" 1 4G 04:00:00 "" "$(r_cmd 01_arm1_inputs.R)")
    record 1 01_arm1_inputs.R "$IN"
    LD=$(sbatch_step greml-ldscore "" 1 16G 06:00:00 \
        "$(task_ranges "${RUN_DIR}/config/ldscore-tasks.tsv" 1)" "$(r_cmd 02_arm1_ldscore.R)")
    record 2 02_arm1_ldscore.R "$LD"
    STRAT_DEPS=""
    REML_JOBS=""
    for CLS in small large; do
        S_RANGE=$(task_ranges "${RUN_DIR}/config/stratify-tasks.tsv" "$CLS")
        [ -z "$S_RANGE" ] && continue
        if [ "$CLS" = small ]; then S_RES=(2 16G 04:00:00); R_RES=(1 8G 08:00:00)
        else S_RES=(4 48G 12:00:00); R_RES=(2 48G 24:00:00); fi
        S=$(sbatch_step "greml-grm-$CLS" "afterok:${LD}" "${S_RES[@]}" "$S_RANGE" \
            "$(r_cmd 03_arm1_stratify_grm.R)")
        record 3 "03_arm1_stratify_grm.R ($CLS)" "$S"
        STRAT_DEPS="${STRAT_DEPS}:${S}"
        R_RANGE=$(task_ranges "${RUN_DIR}/config/reml-tasks.tsv" "$CLS")
        R=$(sbatch_step "greml-reml-$CLS" "afterok:${S}" "${R_RES[@]}" "$R_RANGE" \
            "$(r_cmd 04_arm1_reml.R)")
        record 4 "04_arm1_reml.R ($CLS)" "$R"
        REML_JOBS="${REML_JOBS}:${R}"
    done
    DEPS="afterok:${IN}${REML_JOBS}"
else
    R=$(sbatch_step greml-reml "" 1 8G 06:00:00 \
        "$(task_ranges "${RUN_DIR}/config/locus-tasks.tsv" 1)" "$(r_cmd 05_arm2_greml.R)")
    record 5 05_arm2_greml.R "$R"
    DEPS="afterok:${R}"
fi

SUM=$(sbatch_step greml-summary "$DEPS" 1 16G 02:00:00 "" "$(r_cmd 06_summarize.R)")
record 6 06_summarize.R "$SUM"
GATE=$(sbatch_step greml-gate "afterok:${SUM}" 1 4G 00:30:00 "" "$(r_cmd 07_apply_gate.R)")
record 7 07_apply_gate.R "$GATE"
PLOT=$(sbatch_step greml-plot "afterok:${GATE}" 1 8G 00:30:00 "" "$(r_cmd 08_plot.R)")
record 8 08_plot.R "$PLOT"
FIN=$(sbatch_step greml-final "afterok:${PLOT}" 1 4G 01:00:00 "" "$(r_cmd 09_finalize_run.R)")
record 9 09_finalize_run.R "$FIN"

log_message "job graph recorded in ${JOBS_TSV}"
cat "$JOBS_TSV"
