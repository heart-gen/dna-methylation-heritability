#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
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
H_DIR="${REPO_DIR}/04_repeat_repressive_architecture/_h"

log_job_info
run_r "${H_DIR}/02_test_association.R" --run-id "$RUN_ID"
