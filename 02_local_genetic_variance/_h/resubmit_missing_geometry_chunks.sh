#!/bin/bash

## Resubmit the geometry chunks whose task rows are missing.
##
## Usage: resubmit_missing_geometry_chunks.sh GEOMETRY_RUN_DIR [COMBINE_JOB_ID]
##
## The locus-geometry array writes one row file per VMR, so a chunk that was
## killed leaves exactly its 25 files absent. This script reads the run's own
## task and chunk manifests, finds every task_id with no row file, maps those
## back to their chunk_ids and resubmits only those chunks -- with the faulty
## node excluded, since the recurring cause is a node that root-cancels tasks
## seconds after they start (see submit_observed_regime_cell.sh).
##
## Pass COMBINE_JOB_ID to have the held step_10 job's dependency extended to
## the resubmitted array and released, so the combine cannot run before the
## replacement rows exist. Without it the job id is printed and the caller
## wires the dependency by hand.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 GEOMETRY_RUN_DIR [COMBINE_JOB_ID]" >&2
    exit 2
fi
RUN_DIR=$(cd "$1" && pwd)
COMBINE_JOB=${2:-}

H_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd "${H_DIR}/../.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
PARTITION=${LGV_PARTITION:-short}
MAX_CONCURRENT=${LGV_MAX_CONCURRENT:-200}
EXCLUDE_NODES=${LGV_EXCLUDE_NODES-qnode0287}
EXCLUDE_OPT=()
if [[ -n "${EXCLUDE_NODES}" ]]; then
    EXCLUDE_OPT=(--exclude="${EXCLUDE_NODES}")
fi

TASKS="${RUN_DIR}/config/task-manifest.tsv"
CHUNKS="${RUN_DIR}/config/chunk-manifest.tsv"
ROWS="${RUN_DIR}/results/task_rows"
for path in "${TASKS}" "${CHUNKS}"; do
    [[ -f "${path}" ]] || { echo "Missing manifest: ${path}" >&2; exit 1; }
done

## chunk-manifest.tsv is chunk_id <tab> task_id, one row per task.
MISSING_CHUNKS=$(awk -F '\t' -v rows="${ROWS}" '
    NR == 1 {
        for (i = 1; i <= NF; i++) { col[$i] = i }
        if (!("chunk_id" in col) || !("task_id" in col)) {
            print "chunk-manifest.tsv lacks chunk_id/task_id" > "/dev/stderr"
            exit 1
        }
        next
    }
    {
        f = sprintf("%s/vmr-%07d.tsv", rows, $(col["task_id"]))
        ## <= 0 rather than < 0: an existing but empty row file is as
        ## missing as an absent one, and a real row file has a header.
        if ((getline line < f) <= 0) { bad[$(col["chunk_id"])] = 1 }
        close(f)
    }
    END { for (c in bad) print c }
' "${CHUNKS}" | sort -n | uniq | paste -sd, -)

if [[ -z "${MISSING_CHUNKS}" ]]; then
    echo "No missing task rows in ${RUN_DIR}"
    if [[ -n "${COMBINE_JOB}" ]]; then
        scontrol release "${COMBINE_JOB}"
        echo "Released combine job ${COMBINE_JOB} unchanged"
    fi
    exit 0
fi

COHORT=$(awk -F '\t' '$1 == "cohort" {print $2}' "${RUN_DIR}/manifest.tsv")
REGION=$(awk -F '\t' '$1 == "region" {print $2}' "${RUN_DIR}/manifest.tsv")
EXPORTS="ALL,LGV_GEOMETRY_DIR=${RUN_DIR},LGV_H_DIR=${H_DIR}"
EXPORTS="${EXPORTS},CAL_H2_ENV=${ENV_PATH},V2_REPO_ROOT=${REPO_DIR}"

JOB=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --array="${MISSING_CHUNKS}%${MAX_CONCURRENT}" \
    --job-name="lgv_${COHORT}_${REGION}_geometry_rerun" \
    --output="${RUN_DIR}/logs/%x.%A_%a.log" \
    --export="${EXPORTS}" "${EXCLUDE_OPT[@]}" \
    "${H_DIR}/step_09_locus_geometry.sh" | cut -d';' -f1)

printf 'geometry_rerun\tstep_09_locus_geometry.sh\t%s\t%s\n' \
    "${JOB}" "${MISSING_CHUNKS}" >> "${RUN_DIR}/submitted-jobs.tsv"

echo "Resubmitted chunks ${MISSING_CHUNKS} as ${JOB}"

if [[ -n "${COMBINE_JOB}" ]]; then
    ## afterany, matching the launcher: a cancelled rerun must still be
    ## reconciled rather than leaving the combine unsatisfiable.
    scontrol update JobId="${COMBINE_JOB}" Dependency="afterany:${JOB}"
    scontrol release "${COMBINE_JOB}"
    echo "Combine job ${COMBINE_JOB} now depends on ${JOB} and is released"
fi
