#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomicsguest
#SBATCH --job-name=lsp-oof
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# A chunk is VMRS_PER_ARRAY_TASK sequential VMRs, the same design 02 adopted
# after single-VMR arrays of 60k tasks were mass-cancelled by the scheduler.
# Nested CV is far heavier per VMR than 02's single fit -- outer folds x repeats
# x alpha grid x inner folds -- so the chunk size is smaller and the wall time
# larger. Re-time the first chunks before trusting these numbers at scale.
#SBATCH --time=04:00:00

# SLURM copies the batch script to /var/spool, so ${BASH_SOURCE[0]} does NOT
# resolve to _h/ at run time. V2_REPO_ROOT is exported by the submit driver;
# fall back to the submit directory for a hand-run sbatch.
V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID=${LSP_RUN_ID:?LSP_RUN_ID must be set}
CHUNK_MANIFEST=${LSP_CHUNK_MANIFEST:?LSP_CHUNK_MANIFEST must be set}
ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set}

log_job_info
require_file "$CHUNK_MANIFEST"

# Each task selects its own rows from the manifest rather than being handed them
# on the command line, so a resubmission of a subset of chunks reproduces the
# identical VMR-to-chunk mapping.
mapfile -t TASK_IDS < <(awk -F '\t' -v chunk="${ARRAY_INDEX}" \
    'NR > 1 && $1 == chunk {print $2}' "$CHUNK_MANIFEST")

if [ ${#TASK_IDS[@]} -eq 0 ]; then
    echo "ERROR: chunk ${ARRAY_INDEX} has no tasks in ${CHUNK_MANIFEST}" >&2
    exit 1
fi

log_message "chunk ${ARRAY_INDEX}: ${#TASK_IDS[@]} VMR task(s)"

IFS=,; TASK_LIST="${TASK_IDS[*]}"; unset IFS

run_r "${REPO_DIR}/03_local_snp_prediction/_h/02_fit_oof.R" \
    --run-id "$RUN_ID" \
    --task-ids "$TASK_LIST" \
    ${LSP_EXTRA_ARGS:-}

log_message "chunk ${ARRAY_INDEX} complete"
