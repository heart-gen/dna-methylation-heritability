#!/bin/bash
#
# Resume one or more incomplete observed runs without restoring the original
# one-VMR-per-array-task launch pattern. Each recovery array task receives an
# explicit chunk of pending VMR task IDs and processes them sequentially.
#
# Usage:
#   ../_h/resubmit_pending.sh OBSERVED_RUN_ID [OBSERVED_RUN_ID...]
#
# Environment:
#   VMRS_PER_ARRAY_TASK=25  pending VMRs processed sequentially per array task
#   MAX_CONCURRENT=50       array throttle per cell
#   RESUME_ID=resume-...    optional unique provenance label
#   DRY_RUN=1               report pending/chunk counts without writing/submitting

set -euo pipefail

(( $# >= 1 )) || {
    echo "Usage: $0 OBSERVED_RUN_ID [OBSERVED_RUN_ID...]" >&2
    exit 1
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
VMRS_PER_ARRAY_TASK=${VMRS_PER_ARRAY_TASK:-25}
MAX_CONCURRENT=${MAX_CONCURRENT:-50}
RESUME_ID=${RESUME_ID:-resume-$(date -u +%Y%m%dT%H%M%SZ)}

for pair in "VMRS_PER_ARRAY_TASK:${VMRS_PER_ARRAY_TASK}" \
            "MAX_CONCURRENT:${MAX_CONCURRENT}"; do
    name=${pair%%:*}
    value=${pair#*:}
    if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
        echo "${name} must be a positive integer, got: ${value}" >&2
        exit 1
    fi
done
if [[ ! "${RESUME_ID}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "RESUME_ID contains unsupported characters: ${RESUME_ID}" >&2
    exit 1
fi

pending_for() {
    local root=$1 n=$2 done_file=$3
    find "${root}/results" \
        \( -path '*/summary/vmr-*.tsv' -o -path '*/qc_failures/vmr-*.tsv' \
           -o -path '*/excluded/vmr-*.tsv' \) -type f 2>/dev/null \
        | sed 's#.*/vmr-0*\([0-9]*\)\.tsv#\1#' | sort -n > "${done_file}.all"
    if [[ -s "${done_file}.all" ]] && uniq -d "${done_file}.all" | grep -q .; then
        echo "A VMR appears in multiple non-retryable terminal categories under ${root}" >&2
        return 2
    fi
    uniq "${done_file}.all" > "${done_file}"
    awk 'NR == FNR {done[$1] = 1; next} !($1 in done)' \
        "${done_file}" <(seq 1 "${n}")
}

submitted=0
for run_id in "$@"; do
    ROOT=${ANALYSIS_DIR}/_m/runs/${run_id}
    [[ -d "${ROOT}" ]] || { echo "No such observed run: ${ROOT}" >&2; exit 1; }
    EXPECTED=${ROOT}/config/expected-tasks.tsv
    METADATA=${ROOT}/provenance/run-metadata.tsv
    [[ -s "${EXPECTED}" && -s "${METADATA}" ]] || {
        echo "Run lacks expected-task or metadata manifest: ${ROOT}" >&2
        exit 1
    }
    read -r REGION COHORT N_TASKS < <(awk -F'\t' 'NR == 2 {print $1, $2, $3}' "${EXPECTED}")
    VMR_RUN_DIR=$(awk -F'\t' '$1 == "upstream_vmr_run_dir" {print $2}' "${METADATA}")
    RUN_SCRIPT_DIR=${ROOT}/code/_h
    RUN_MODEL=${ROOT}/config/elastic-net-calibration.rds
    [[ -s "${RUN_SCRIPT_DIR}/step_5_estimate_observed_vmr.sh" && -s "${RUN_MODEL}" ]] || {
        echo "Run code snapshot or calibration model is missing: ${ROOT}" >&2
        exit 1
    }

    done_file=$(mktemp "${TMPDIR:-/tmp}/cal-h2-done.XXXXXX")
    pending_file=$(mktemp "${TMPDIR:-/tmp}/cal-h2-pending.XXXXXX")
    trap 'rm -f "${done_file:-}" "${done_file:-}.all" "${pending_file:-}"' EXIT
    pending_for "${ROOT}" "${N_TASKS}" "${done_file}" > "${pending_file}"
    N_PENDING=$(wc -l < "${pending_file}")
    if (( N_PENDING == 0 )); then
        printf '%-46s complete; no pending VMRs\n' "${run_id}"
        rm -f "${done_file}" "${done_file}.all" "${pending_file}"
        trap - EXIT
        continue
    fi
    N_CHUNKS=$(( (N_PENDING + VMRS_PER_ARRAY_TASK - 1) / VMRS_PER_ARRAY_TASK ))
    MANIFEST=${ROOT}/config/${RESUME_ID}-chunk-manifest.tsv
    if [[ -e "${MANIFEST}" ]]; then
        echo "Resume manifest already exists; choose a new RESUME_ID: ${MANIFEST}" >&2
        exit 1
    fi

    if [[ "${DRY_RUN:-0}" == 1 ]]; then
        printf '%-46s %6d pending VMRs -> %d chunks (size <= %d, throttle %d)\n' \
            "${run_id}" "${N_PENDING}" "${N_CHUNKS}" \
            "${VMRS_PER_ARRAY_TASK}" "${MAX_CONCURRENT}"
        rm -f "${done_file}" "${done_file}.all" "${pending_file}"
        trap - EXIT
        continue
    fi

    awk -v size="${VMRS_PER_ARRAY_TASK}" 'BEGIN {OFS="\t"; print "chunk_id", "task_id"}
        {print int((NR - 1) / size) + 1, $1}' "${pending_file}" > "${MANIFEST}"
    JOB=$(sbatch --parsable --account="${ACCOUNT}" \
        --array="1-${N_CHUNKS}%${MAX_CONCURRENT}" \
        --job-name="cal_h2_${COHORT}_${REGION}" \
        --output="${ROOT}/logs/%x.%A_%a.log" \
        --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_REPO_ROOT="${REPO_ROOT}",CAL_H2_CALIBRATION_MODEL="${RUN_MODEL}",CAL_H2_OBSERVED_OUTPUT_ROOT="${ROOT}/results",CAL_H2_VMR_RUN_DIR="${VMR_RUN_DIR}",CAL_H2_COHORT="${COHORT}",CAL_H2_CHUNK_MANIFEST="${MANIFEST}",CAL_H2_CHUNK_SET_ID="${RESUME_ID}",CAL_H2_PLINK_ROOT=,CAL_H2_PHENOTYPE_ROOT=,CAL_H2_RECOVERED_PLINK_ROOT=,CAL_H2_WRITE_DIAGNOSTICS=FALSE,REGION="${REGION}",POPULATION="${COHORT}" \
        "${RUN_SCRIPT_DIR}/step_5_estimate_observed_vmr.sh")
    printf 'resume:%s\t%s\t%s\t%s\t%s\n' "${RESUME_ID}" "${REGION}" "${JOB}" \
        "${N_CHUNKS}" "${N_PENDING}" >> "${ROOT}/provenance/submitted-jobs.tsv"
    printf '%-46s %6d pending VMRs -> job %s (%d chunks, size <= %d)\n' \
        "${run_id}" "${N_PENDING}" "${JOB}" "${N_CHUNKS}" "${VMRS_PER_ARRAY_TASK}"
    submitted=$((submitted + 1))

    rm -f "${done_file}" "${done_file}.all" "${pending_file}"
    trap - EXIT
done

if [[ "${DRY_RUN:-0}" == 1 ]]; then
    echo "Dry run only; no manifests were written and no jobs were submitted."
elif (( submitted == 0 )); then
    echo "All requested runs are already complete."
fi
