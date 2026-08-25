#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=lsp_combine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.err
#
# Module 03 step 3: reconcile every task and pool the out-of-fold predictions.
#
# Env-driven with no positional arguments, so any single stage can be
# resubmitted by hand against an existing run without editing this file:
#   LSP_RUN_ID=<id> sbatch _h/step_3_combine_oof.sh

# SLURM copies the batch script to /var/spool, so ${BASH_SOURCE[0]} does NOT
# resolve to _h/ at run time. V2_REPO_ROOT is exported by the submit driver;
# fall back to the submit directory for a hand-run sbatch.
V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID=${LSP_RUN_ID:?LSP_RUN_ID must be set}
H_DIR="${REPO_DIR}/03_local_snp_prediction/_h"
EXTRA=${LSP_EXTRA_ARGS:-}

log_job_info
run_r "${H_DIR}/03_combine_oof.R" --run-id "$RUN_ID" $EXTRA
