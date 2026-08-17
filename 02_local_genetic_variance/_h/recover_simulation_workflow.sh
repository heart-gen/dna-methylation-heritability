#!/bin/bash
# One bounded reconciliation/recovery after all original simulation arrays are
# terminal. It submits exactly missing scenario IDs, then a new fit and
# independent-evaluation chain. It never loops or resubmits itself.

set -euo pipefail
if (( $# != 1 )); then
    echo "Usage: $0 CALIBRATION_RUN_ID" >&2
    exit 1
fi
RUN_ID=$1
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
RUN_ROOT=${ANALYSIS_DIR}/_m/runs/${RUN_ID}
MANIFEST=${RUN_ROOT}/config/scenarios.tsv
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
MAX_CONCURRENT=${MAX_CONCURRENT:-20}
SCENARIOS_PER_ARRAY_TASK=${SCENARIOS_PER_ARRAY_TASK:-10}
if [[ ! "${SCENARIOS_PER_ARRAY_TASK}" =~ ^[1-9][0-9]*$ ]]; then
    echo "SCENARIOS_PER_ARRAY_TASK must be a positive integer" >&2
    exit 1
fi
[[ -s "${MANIFEST}" ]] || { echo "Missing scenario manifest: ${MANIFEST}" >&2; exit 1; }

original_jobs=$(awk -F '\t' '$1 == "simulation" {print $2}' \
    "${RUN_ROOT}/provenance/submitted-jobs.tsv" | paste -sd, -)
if [[ -n "${original_jobs}" ]] && [[ -n "$(squeue -h -j "${original_jobs}")" ]]; then
    echo "Original simulation arrays are still active; recovery is premature" >&2
    exit 1
fi

RECOVERY_TAG=$(date -u +%Y%m%dT%H%M%SZ)
RECOVERY_DIR=${RUN_ROOT}/recovery/${RECOVERY_TAG}
RECOVERY_SCRIPT_DIR=${RUN_ROOT}/code/recovery-${RECOVERY_TAG}/_h
mkdir -p "${RECOVERY_DIR}" "$(dirname "${RECOVERY_SCRIPT_DIR}")"
cp -a "${SCRIPT_DIR}" "${RECOVERY_SCRIPT_DIR}"

EXPECTED=$(mktemp "${TMPDIR:-/tmp}/cal-h2-expected.XXXXXX")
COMPLETED=$(mktemp "${TMPDIR:-/tmp}/cal-h2-completed.XXXXXX")
MISSING=$(mktemp "${TMPDIR:-/tmp}/cal-h2-missing.XXXXXX")
trap 'rm -f "${EXPECTED}" "${COMPLETED}" "${MISSING}"' EXIT
awk -F '\t' 'NR > 1 {print $1}' "${MANIFEST}" | sort -n > "${EXPECTED}"
find "${RUN_ROOT}/raw" -type f -name 'scenario-*.tsv' -printf '%f\n' \
    | sed -E 's/scenario-0*([0-9]+)\.tsv/\1/' | sort -n -u > "${COMPLETED}"
comm -23 "${EXPECTED}" "${COMPLETED}" > "${MISSING}"

RECOVERY_MANIFEST=${RECOVERY_DIR}/missing-scenarios.tsv
{
    printf 'chunk_id\tscenario_id\n'
    awk -v size="${SCENARIOS_PER_ARRAY_TASK}" \
        '{print int((NR - 1) / size) + 1 "\t" $1}' "${MISSING}"
} > "${RECOVERY_MANIFEST}"
MISSING_N=$(wc -l < "${MISSING}")
RECOVERY_CHUNKS=$(( (MISSING_N + SCENARIOS_PER_ARRAY_TASK - 1) / SCENARIOS_PER_ARRAY_TASK ))
EXPECTED_N=$(wc -l < "${EXPECTED}")
COMPLETED_N=$(wc -l < "${COMPLETED}")
printf 'expected_scenarios\tcompleted_before_recovery\tmissing_scenarios\n%s\t%s\t%s\n' \
    "${EXPECTED_N}" "${COMPLETED_N}" "${MISSING_N}" > \
    "${RECOVERY_DIR}/reconciliation-before-recovery.tsv"

RECOVERY_JOB_IDS=()
if (( MISSING_N > 0 )); then
    MAX_ARRAY_TASKS=9000
    for (( OFFSET=0; OFFSET<RECOVERY_CHUNKS; OFFSET+=MAX_ARRAY_TASKS )); do
        REMAINING=$((RECOVERY_CHUNKS - OFFSET))
        BATCH_TASKS=$((REMAINING < MAX_ARRAY_TASKS ? REMAINING : MAX_ARRAY_TASKS))
        START=$((OFFSET + 1)); END=$((OFFSET + BATCH_TASKS))
        job=$(sbatch --parsable --account="${ACCOUNT}" \
            --array="${START}-${END}%${MAX_CONCURRENT}" \
            --output="${RUN_ROOT}/logs/%x.%A_%a.log" \
            --export=ALL,CAL_H2_SCRIPT_DIR="${RECOVERY_SCRIPT_DIR}",SCENARIO_MANIFEST="${MANIFEST}",CAL_H2_RECOVERY_MANIFEST="${RECOVERY_MANIFEST}",SIMULATION_OUTPUT_ROOT="${RUN_ROOT}/raw" \
            "${RECOVERY_SCRIPT_DIR}/step_2_recover_simulation.sh")
        RECOVERY_JOB_IDS+=("${job%%;*}")
    done
    RECOVERY_DEPENDENCY=$(IFS=:; echo "${RECOVERY_JOB_IDS[*]}")
else
    RECOVERY_DEPENDENCY=""
fi

dependency_args=()
if [[ -n "${RECOVERY_DEPENDENCY}" ]]; then
    dependency_args+=(--dependency="afterok:${RECOVERY_DEPENDENCY}")
fi
CAL_JOB=$(sbatch --parsable --account="${ACCOUNT}" "${dependency_args[@]}" \
    --output="${RUN_ROOT}/logs/%x.%j.log" \
    --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RECOVERY_SCRIPT_DIR}",CAL_H2_RUN_ROOT="${RUN_ROOT}" \
    "${RECOVERY_SCRIPT_DIR}/step_3_fit_calibration.sh")
CAL_JOB_ID=${CAL_JOB%%;*}
EVAL_JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --dependency="afterok:${CAL_JOB_ID}" \
    --output="${RUN_ROOT}/logs/%x.%j.log" \
    --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RECOVERY_SCRIPT_DIR}",CAL_H2_RUN_ROOT="${RUN_ROOT}",CAL_H2_FAIL_ON_REJECTION=TRUE \
    "${RECOVERY_SCRIPT_DIR}/step_4_evaluate_calibration.sh")
EVAL_JOB_ID=${EVAL_JOB%%;*}

{
    printf 'stage\tjob_id\tdependency\n'
    for job_id in "${RECOVERY_JOB_IDS[@]}"; do
        printf 'simulation_recovery\t%s\tNA\n' "${job_id}"
    done
    printf 'calibration_recovery\t%s\tafterok:%s\n' "${CAL_JOB_ID}" "${RECOVERY_DEPENDENCY:-NA}"
    printf 'evaluation_recovery\t%s\tafterok:%s\n' "${EVAL_JOB_ID}" "${CAL_JOB_ID}"
} > "${RECOVERY_DIR}/submitted-jobs.tsv"
echo "Missing scenarios: ${MISSING_N} in ${RECOVERY_CHUNKS} chunks; recovery arrays: ${RECOVERY_JOB_IDS[*]:-none}"
echo "Calibration: ${CAL_JOB_ID}; evaluation gate: ${EVAL_JOB_ID}"
