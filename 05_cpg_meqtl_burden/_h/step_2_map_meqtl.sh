#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomicsguest
#SBATCH --job-name=cmb-meqtl
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
# NOTE: no --array here. The driver passes the array spec, so a smoke run can
# map one chromosome without editing this file, while the production default
# (1-22) stays the config autosome policy.
#SBATCH --time=06:00:00
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# One array task per autosome. meqtl_parameters.yml fixes chromosomes:
# autosomal, so there is no chrX task and none is silently skipped -- the array
# bound is the policy.

# SLURM copies the batch script to /var/spool, so ${BASH_SOURCE[0]} does NOT
# resolve to _h/ at run time. V2_REPO_ROOT is exported by the submit driver;
# fall back to the submit directory for a hand-run sbatch.
V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID=${CMB_RUN_ID:?CMB_RUN_ID must be set}
CHROM=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set}
RUN_DIR="${REPO_DIR}/05_cpg_meqtl_burden/_m/runs/${RUN_ID}"

log_job_info
require_file "${RUN_DIR}/results/tested-cpg-membership.tsv"

# tensorqtl is not in the epigenomics env (it drags a pinned torch build). It
# lives in `genomics`, alongside pgenlib, py_qvalue and plink -- verified
# 2026-08-23. There is no separate `tensorqtl` env on this system.
V2_ENV_TENSORQTL="${V2_ENV_TENSORQTL:-$ENV_PATH/genomics}"
require_file "$V2_ENV_TENSORQTL"

# One array task owns one autosome end to end: preparing that chromosome's
# inputs and then mapping it. Keeping the two together means a resubmitted
# chromosome rebuilds exactly the inputs it maps, with no stale-input window.
log_message "preparing tensorqtl inputs for chr${CHROM}"
conda run --no-capture-output -p "$V2_ENV_TENSORQTL" \
    python "${REPO_DIR}/05_cpg_meqtl_burden/_h/01b_prepare_meqtl_inputs.py" \
    --run-id "$RUN_ID" \
    --chrom "$CHROM" \
    --threads "$V2_THREADS"

log_message "mapping cis-meQTL for chr${CHROM}"
conda run --no-capture-output -p "$V2_ENV_TENSORQTL" \
    python "${REPO_DIR}/05_cpg_meqtl_burden/_h/02_map_cpg_meqtl.py" \
    --run-id "$RUN_ID" \
    --chrom "$CHROM" \
    --threads "$V2_THREADS"

log_message "chr${CHROM} complete"
