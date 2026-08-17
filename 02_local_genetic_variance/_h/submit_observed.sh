#!/bin/bash
#
# 02_local_genetic_variance: estimate observed local genetic variance against an
# ACCEPTED 01_vmr_catalog run.
#
# Why this exists rather than submit_observed_calibrated_workflow.sh: that
# script was migrated from calibrated-simulation-analysis/ unchanged, and while
# 04_estimate_observed_vmr.R and step_5_estimate_observed_vmr.sh were re-pointed
# at corrected catalogs, the submitter was not. It still
#
#   - sizes its arrays from vmr-analysis/all_individuals/<region>/_m/vmr.bed,
#     the V1-invalidated legacy catalog;
#   - defaults CAL_H2_PLINK_ROOT to the lead author's tree;
#   - expects the calibration model under _m/runs/<CALIBRATION_RUN_ID>/, which
#     does not exist in this module (the model is frozen under
#     _m/calibration_frozen/).
#
# It is left in place for provenance. This is the active entry point.
#
# Usage, from the module's _m/ directory:
#   cd 02_local_genetic_variance/_m
#   ../_h/submit_observed.sh <OBSERVED_RUN_ID> <AA|all_individuals> <region>
#
# Environment:
#   SMOKE_N=50        estimate only the first N VMRs. Not production.
#   VMRS_PER_ARRAY_TASK  VMRs processed sequentially in one task (default 25).
#   MAX_CONCURRENT       array throttle per cell (default 50).
#   DRY_RUN=1         print the plan without submitting.

set -euo pipefail

usage() {
    echo "Usage: $0 OBSERVED_RUN_ID <AA|all_individuals> <caudate|dlpfc|hippocampus>" >&2
}
if (( $# != 3 )); then usage; exit 1; fi

OBSERVED_RUN_ID=$1
COHORT=$2
REGION=$3

case "$COHORT" in AA|all_individuals) ;; *) usage; exit 1 ;; esac
case "$REGION" in caudate|dlpfc|hippocampus) ;; *) usage; exit 1 ;; esac
if [[ ! "${OBSERVED_RUN_ID}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "Run ID contains unsupported characters: ${OBSERVED_RUN_ID}" >&2
    exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
VMRS_PER_ARRAY_TASK=${VMRS_PER_ARRAY_TASK:-25}
MAX_CONCURRENT=${MAX_CONCURRENT:-50}
for pair in "VMRS_PER_ARRAY_TASK:${VMRS_PER_ARRAY_TASK}" \
            "MAX_CONCURRENT:${MAX_CONCURRENT}"; do
    name=${pair%%:*}
    value=${pair#*:}
    if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
        echo "${name} must be a positive integer, got: ${value}" >&2
        exit 1
    fi
done

# The accepted upstream catalog for this cell. The gate that made it accepted is
# recorded in 01_vmr_catalog/README.md (AGENTS.md 6); this script asserts the
# run is sealed and complete, not that the gate was reviewed.
VMR_RUN_DIR=${CAL_H2_VMR_RUN_DIR:-${REPO_ROOT}/01_vmr_catalog/_m/runs/vmrcat-${COHORT}-${REGION}-20260816}
CALIBRATION_ROOT=${CAL_H2_CALIBRATION_ROOT:-${ANALYSIS_DIR}/_m/calibration_frozen}
if [[ -s "${CALIBRATION_ROOT}/calibration/elastic-net-calibration.rds" ]]; then
    CALIBRATION_MODEL=${CALIBRATION_ROOT}/calibration/elastic-net-calibration.rds
    CALIBRATION_PERFORMANCE=${CALIBRATION_ROOT}/evaluation/calibration-performance-overall.tsv
    CALIBRATION_ACCEPTANCE=${CALIBRATION_ROOT}/evaluation/acceptance-results.tsv
    CALIBRATION_VALIDATION_METADATA=${CALIBRATION_ROOT}/provenance/run-metadata.tsv
else
    CALIBRATION_MODEL=${CALIBRATION_ROOT}/elastic-net-calibration.rds
    CALIBRATION_PERFORMANCE=${CALIBRATION_ROOT}/calibration-performance-overall.tsv
    CALIBRATION_ACCEPTANCE=${CALIBRATION_ROOT}/acceptance-results.tsv
    CALIBRATION_VALIDATION_METADATA=${CALIBRATION_ROOT}/validation-metadata.tsv
fi
OBSERVED_ROOT=${ANALYSIS_DIR}/_m/runs/${OBSERVED_RUN_ID}

for f in "${VMR_RUN_DIR}/manifest.tsv" "${VMR_RUN_DIR}/vmr/vmr.bed" \
         "${VMR_RUN_DIR}/task_reconciliation.tsv" "${CALIBRATION_MODEL}"; do
    [[ -s "$f" ]] || { echo "Required input is missing: $f" >&2; exit 1; }
done
[[ -d "${VMR_RUN_DIR}/plink_format" ]] || {
    echo "Upstream run has no plink_format/; run 01_vmr_catalog step_4 first" >&2; exit 1; }

# AGENTS.md 14: a run submitted with --allow-unlocked is never citable as
# production, and a non-zero failed count in the upstream reconciliation is a
# blocking dependency gate (AGENTS.md 6). Check both rather than trusting the
# README.
man_field() { awk -v f="$1" -F'\t' '$1 == f {print $2}' "${VMR_RUN_DIR}/manifest.tsv"; }
UPSTREAM_SMOKE=$(man_field smoke_run)
UPSTREAM_VMR_SET_ID=$(man_field vmr_set_id)
if [[ "${UPSTREAM_SMOKE}" == "TRUE" ]]; then
    echo "Upstream ${VMR_RUN_DIR} is a smoke run; not citable as production" >&2
    exit 1
fi
[[ -n "${UPSTREAM_VMR_SET_ID}" ]] || {
    echo "Upstream run is not sealed: manifest has no vmr_set_id" >&2; exit 1; }
UPSTREAM_FAILED=$(awk -F'\t' '$1 ~ /^(failed|unaccounted|unexpected|qc_failed)$/ {s += $2} END {print s+0}' \
    "${VMR_RUN_DIR}/task_reconciliation.tsv")
if (( UPSTREAM_FAILED != 0 )); then
    echo "Upstream reconciliation reports ${UPSTREAM_FAILED} unexplained tasks; blocked" >&2
    exit 1
fi

for f in "${CALIBRATION_PERFORMANCE}" "${CALIBRATION_ACCEPTANCE}"; do
    [[ -s "${f}" ]] || { echo "Calibration gate input is missing: ${f}" >&2; exit 1; }
done
ACTUAL_MODEL_SHA=$(sha256sum "${CALIBRATION_MODEL}" | awk '{print $1}')

# Fail-closed acceptance gate. The estimator may not be applied to observed data
# unless every criterion in config/acceptance-criteria.tsv passes.
# write_tsv() writes to a temp file and renames, so --output cannot be /dev/null.
GATE_FILE=$(mktemp "${TMPDIR:-/tmp}/cal-h2-gate.XXXXXX.tsv")
trap 'rm -f "${GATE_FILE}"' EXIT
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/06_check_acceptance.R" \
    --performance="${CALIBRATION_PERFORMANCE}" \
    --criteria="${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
    --model="${CALIBRATION_MODEL}" \
    --output="${GATE_FILE}" \
    --fail-on-rejection=TRUE

N_VMRS=$(wc -l < "${VMR_RUN_DIR}/vmr/vmr.bed")
N_TASKS=${N_VMRS}
if [[ -n "${SMOKE_N:-}" ]]; then
    if [[ ! "${SMOKE_N}" =~ ^[1-9][0-9]*$ ]]; then
        echo "SMOKE_N must be a positive integer, got: ${SMOKE_N}" >&2
        exit 1
    fi
    N_TASKS=$(( SMOKE_N < N_VMRS ? SMOKE_N : N_VMRS ))
fi
N_ARRAY_TASKS=$(( (N_TASKS + VMRS_PER_ARRAY_TASK - 1) / VMRS_PER_ARRAY_TASK ))

if [[ "${DRY_RUN:-0}" == "1" ]]; then
    cat >&2 <<EOF
[dry-run] observed run : ${OBSERVED_ROOT}
[dry-run] upstream     : ${VMR_RUN_DIR}
[dry-run] vmr_set_id   : ${UPSTREAM_VMR_SET_ID}
[dry-run] VMRs         : ${N_TASKS} of ${N_VMRS}
[dry-run] array        : 1-${N_ARRAY_TASKS}%${MAX_CONCURRENT}
[dry-run] chunk size   : ${VMRS_PER_ARRAY_TASK} sequential VMRs/task
EOF
    exit 0
fi

if [[ -e "${OBSERVED_ROOT}" ]]; then
    echo "Observed run already exists: ${OBSERVED_ROOT}" >&2
    exit 1
fi
mkdir -p "${OBSERVED_ROOT}"/{config,code,logs,provenance,results/chunk_status}
cp -a "${SCRIPT_DIR}" "${OBSERVED_ROOT}/code/_h"
RUN_SCRIPT_DIR=${OBSERVED_ROOT}/code/_h
cp "${CALIBRATION_MODEL}" "${OBSERVED_ROOT}/config/elastic-net-calibration.rds"
cp "${ANALYSIS_DIR}/config/acceptance-criteria.tsv" \
   "${CALIBRATION_ACCEPTANCE}" "${CALIBRATION_PERFORMANCE}" \
   "${OBSERVED_ROOT}/config/"
if [[ -s "${CALIBRATION_VALIDATION_METADATA}" ]]; then
    cp "${CALIBRATION_VALIDATION_METADATA}" \
       "${OBSERVED_ROOT}/config/calibration-validation-metadata.tsv"
fi
cp "${GATE_FILE}" "${OBSERVED_ROOT}/config/calibration-acceptance-results.tsv"
RUN_MODEL=${OBSERVED_ROOT}/config/elastic-net-calibration.rds

# One row per VMR makes the assignment auditable and supports recovery manifests
# containing arbitrary pending task IDs. Array tasks select all rows matching
# their chunk_id and process those VMRs sequentially.
CHUNK_MANIFEST=${OBSERVED_ROOT}/config/chunk-manifest.tsv
awk -v n="${N_TASKS}" -v size="${VMRS_PER_ARRAY_TASK}" 'BEGIN {
    OFS="\t"; print "chunk_id", "task_id"
    for (task = 1; task <= n; task++) print int((task - 1) / size) + 1, task
}' > "${CHUNK_MANIFEST}"
awk -F'\t' -v n="${N_TASKS}" -v chunks="${N_ARRAY_TASKS}" '
    NR == 1 {if ($1 != "chunk_id" || $2 != "task_id") exit 2; next}
    {rows++; seen[$2]++; if ($1 < 1 || $1 > chunks) exit 3}
    END {
        if (rows != n) exit 4
        for (i = 1; i <= n; i++) if (seen[i] != 1) exit 5
    }
' "${CHUNK_MANIFEST}" || {
    echo "Internal error: generated chunk manifest failed validation" >&2
    exit 1
}

printf 'region\tpopulation\texpected_tasks\texpected_chunks\tvmrs_per_array_task\n' > \
    "${OBSERVED_ROOT}/config/expected-tasks.tsv"
printf '%s\t%s\t%s\t%s\t%s\n' "${REGION}" "${COHORT}" "${N_TASKS}" \
    "${N_ARRAY_TASKS}" "${VMRS_PER_ARRAY_TASK}" >> \
    "${OBSERVED_ROOT}/config/expected-tasks.tsv"

{
    printf 'field\tvalue\n'
    printf 'observed_run_id\t%s\n' "${OBSERVED_RUN_ID}"
    printf 'cohort\t%s\n' "${COHORT}"
    printf 'region\t%s\n' "${REGION}"
    printf 'upstream_vmr_run_id\t%s\n' "$(basename "${VMR_RUN_DIR}")"
    printf 'upstream_vmr_set_id\t%s\n' "${UPSTREAM_VMR_SET_ID}"
    printf 'upstream_vmr_run_dir\t%s\n' "${VMR_RUN_DIR}"
    printf 'calibration_model_sha256\t%s\n' "${ACTUAL_MODEL_SHA}"
    printf 'calibration_root\t%s\n' "${CALIBRATION_ROOT}"
    printf 'n_vmrs_in_catalog\t%s\n' "${N_VMRS}"
    printf 'n_vmrs_submitted\t%s\n' "${N_TASKS}"
    printf 'n_array_tasks_submitted\t%s\n' "${N_ARRAY_TASKS}"
    printf 'vmrs_per_array_task\t%s\n' "${VMRS_PER_ARRAY_TASK}"
    printf 'max_concurrent_array_tasks\t%s\n' "${MAX_CONCURRENT}"
    printf 'smoke_run\t%s\n' "$([[ -n "${SMOKE_N:-}" ]] && echo TRUE || echo FALSE)"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "${OBSERVED_ROOT}/provenance/run-metadata.tsv"
git -C "${REPO_ROOT}" rev-parse HEAD > "${OBSERVED_ROOT}/provenance/git-commit.txt"

cd "${REPO_ROOT}"
# CAL_H2_PLINK_ROOT is deliberately EMPTY: a non-empty value makes the adapter
# resolve genotypes as <plink_root>/<region>/_m/plink_format, the legacy layout.
# Empty makes it read plink_format/ inside the accepted catalog run.
JOB=$(sbatch --parsable --account="${ACCOUNT}" \
    --array="1-${N_ARRAY_TASKS}%${MAX_CONCURRENT}" \
    --job-name="cal_h2_${COHORT}_${REGION}" \
    --output="${OBSERVED_ROOT}/logs/%x.%A_%a.log" \
    --export=ALL,CAL_H2_ANALYSIS_DIR="${ANALYSIS_DIR}",CAL_H2_SCRIPT_DIR="${RUN_SCRIPT_DIR}",CAL_H2_REPO_ROOT="${REPO_ROOT}",CAL_H2_CALIBRATION_MODEL="${RUN_MODEL}",CAL_H2_EXPECTED_MODEL_SHA256="${ACTUAL_MODEL_SHA}",CAL_H2_OBSERVED_OUTPUT_ROOT="${OBSERVED_ROOT}/results",CAL_H2_VMR_RUN_DIR="${VMR_RUN_DIR}",CAL_H2_COHORT="${COHORT}",CAL_H2_CHUNK_MANIFEST="${CHUNK_MANIFEST}",CAL_H2_CHUNK_SET_ID=initial,CAL_H2_PLINK_ROOT=,CAL_H2_PHENOTYPE_ROOT=,CAL_H2_RECOVERED_PLINK_ROOT=,CAL_H2_WRITE_DIAGNOSTICS=FALSE,REGION="${REGION}",POPULATION="${COHORT}" \
    "${RUN_SCRIPT_DIR}/step_5_estimate_observed_vmr.sh")

printf 'stage\tregion\tjob_id\tarray_tasks\tvmrs\nestimate\t%s\t%s\t%s\t%s\n' \
    "${REGION}" "${JOB}" "${N_ARRAY_TASKS}" "${N_TASKS}" > \
    "${OBSERVED_ROOT}/provenance/submitted-jobs.tsv"

cat <<EOF

Observed run : ${OBSERVED_ROOT}
Upstream     : $(basename "${VMR_RUN_DIR}")  (${UPSTREAM_VMR_SET_ID})
Array        : ${JOB}  1-${N_ARRAY_TASKS}%${MAX_CONCURRENT}
Work         : ${N_TASKS} of ${N_VMRS} VMRs (${VMRS_PER_ARRAY_TASK} sequential VMRs/task)
EOF
