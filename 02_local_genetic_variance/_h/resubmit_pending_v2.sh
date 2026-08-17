#!/bin/bash
#
# 02_local_genetic_variance: finish an observed run by resubmitting only the
# tasks that have no result yet.
#
# Why this exists: the first full v2 submission put 60,801 array tasks in the
# queue at once (six cells, one array each). 10,801 completed and the remaining
# 50,000 were CANCELLED by uid 0 -- an administrative or scheduler action, not a
# job failure. None of the cancelled tasks started or wrote a log, and every
# task that did run produced output. The cause was not established; the working
# assumption is that a 60k-task burst is more than the site tolerates from one
# user, so this driver keeps the outstanding footprint small and tops it up.
#
# It is idempotent. A task counts as done if it wrote a summary, a QC failure,
# or an exclusion, so rerunning it after a partial round only picks up the gaps
# and never recomputes a locus that already has a result.
#
# Usage, from the module's _m/ directory:
#   ../_h/resubmit_pending_v2.sh lgv-AA-caudate-20260816 [more run ids...]
#
# Environment:
#   CHUNK=2000            tasks submitted per cell per round
#   MAX_CONCURRENT=300    running tasks per chunk
#   POLL=120              seconds between rounds

set -euo pipefail

(( $# >= 1 )) || { echo "Usage: $0 OBSERVED_RUN_ID [OBSERVED_RUN_ID...]" >&2; exit 1; }

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
CHUNK=${CHUNK:-2000}
MAX_CONCURRENT=${MAX_CONCURRENT:-300}
POLL=${POLL:-120}

# A task is done when it has a result of ANY kind. summary, qc_failures and
# excluded are all terminal outcomes; only the absence of all three means the
# locus was never estimated.
pending_for() {
    local root=$1 n=$2
    find "${root}/results" \
        \( -path '*/summary/*.tsv' -o -path '*/qc_failures/*.tsv' \
           -o -path '*/excluded/*.tsv' \) 2>/dev/null \
        | sed 's#.*/vmr-0*\([0-9]*\)\.tsv#\1#' | sort -u > "${root}/config/done-tasks.txt"
    comm -23 <(seq 1 "$n" | sort) <(sort "${root}/config/done-tasks.txt") | sort -n
}

round=0
while :; do
    round=$((round + 1))
    submitted=0
    jobs_this_round=()
    for run_id in "$@"; do
        ROOT=${ANALYSIS_DIR}/_m/observed-runs/${run_id}
        [[ -d "$ROOT" ]] || { echo "No such observed run: $ROOT" >&2; exit 1; }
        read -r REGION COHORT N_TASKS < <(awk -F'\t' 'NR==2{print $1, $2, $3}' \
            "${ROOT}/config/expected-tasks.tsv")
        VMR_RUN_DIR=$(awk -F'\t' '$1=="upstream_vmr_run_dir"{print $2}' \
            "${ROOT}/provenance/run-metadata.tsv")
        RUN_SCRIPT_DIR=${ROOT}/code/_h
        RUN_MODEL=${ROOT}/config/elastic-net-calibration.rds

        mapfile -t pending < <(pending_for "$ROOT" "$N_TASKS")
        (( ${#pending[@]} > 0 )) || continue

        # The manifest is read by step_5 at row (array index + 1), so it needs a
        # header line. One manifest per round, kept for provenance.
        MANIFEST=${ROOT}/config/pending-round${round}.tsv
        printf 'task_id\n' > "$MANIFEST"
        printf '%s\n' "${pending[@]:0:${CHUNK}}" >> "$MANIFEST"
        N_CHUNK=$(( ${#pending[@]} < CHUNK ? ${#pending[@]} : CHUNK ))

        JOB=$(sbatch --parsable --account="${ACCOUNT}" \
            --array="1-${N_CHUNK}%${MAX_CONCURRENT}" \
            --job-name="cal_h2_${COHORT}_${REGION}" \
            --output="${ROOT}/logs/%x.%A_%a.log" \
            --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_REPO_ROOT="${REPO_ROOT}",CAL_H2_CALIBRATION_MODEL="${RUN_MODEL}",CAL_H2_OBSERVED_OUTPUT_ROOT="${ROOT}/results",CAL_H2_VMR_RUN_DIR="${VMR_RUN_DIR}",CAL_H2_COHORT="${COHORT}",CAL_H2_TASK_MANIFEST="${MANIFEST}",CAL_H2_WRITE_DIAGNOSTICS=FALSE,REGION="${REGION}",POPULATION="${COHORT}" \
            "${RUN_SCRIPT_DIR}/step_5_estimate_observed_vmr.sh")
        jobs_this_round+=("$JOB")
        printf 'resubmit_round%s\t%s\t%s\t%s\n' "$round" "$REGION" "$JOB" "$N_CHUNK" \
            >> "${ROOT}/provenance/submitted-jobs.tsv"
        printf '%-46s round %-3s %6d pending  ->  job %s (%d tasks)\n' \
            "$run_id" "$round" "${#pending[@]}" "$JOB" "$N_CHUNK"
        submitted=$((submitted + 1))
    done

    (( submitted > 0 )) || { echo "All runs complete."; break; }

    ids=$(IFS=,; echo "${jobs_this_round[*]}")
    while [[ -n "$(squeue -u "$USER" -h -j "$ids" 2>/dev/null)" ]]; do sleep "$POLL"; done
done
