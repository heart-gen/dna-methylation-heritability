#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=rra_assoc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.err
#
# Module 04 step 2: fit the primary model and every locked sensitivity.
#
#   RRA_RUN_ID=<id> sbatch _h/step_2_associate.sh

V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID=${RRA_RUN_ID:?RRA_RUN_ID must be set}
# Analysis code is executed from the run's immutable snapshot under _m, not
# from the live _h/. V2_RUN_CODE is exported by the submit driver and points at
# {run_dir}/code/_h. Without this, editing _h/ while a run is queued silently
# changes what that run executes, and the snapshot and recorded git_commit
# attest to code that never ran. Falls back to live _h/ for a hand-run sbatch.
H_DIR="${V2_RUN_CODE:-${REPO_DIR}/04_repeat_repressive_architecture/_h}"

log_job_info
run_r "${H_DIR}/02_test_association.R" --run-id "$RUN_ID"
