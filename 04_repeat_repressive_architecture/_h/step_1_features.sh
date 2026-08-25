#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=rra_feat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.err
#
# Module 04 step 1: build the per-VMR outcome and covariate table for one cell.
#
# 64G and four hours because two covariates are genuinely heavy: WGBS coverage
# loads one BSseq object per chromosome, and the cell-composition R2 opens one
# phenotype file per VMR. This is not a job for the submit host.
#
#   RRA_RUN_ID=<id> sbatch _h/step_1_features.sh

V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID=${RRA_RUN_ID:?RRA_RUN_ID must be set}
H_DIR="${REPO_DIR}/04_repeat_repressive_architecture/_h"

log_job_info
run_r "${H_DIR}/01_build_features.R" --run-id "$RUN_ID" ${RRA_EXTRA_ARGS:-}
