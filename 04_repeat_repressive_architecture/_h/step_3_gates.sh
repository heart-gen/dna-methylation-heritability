#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=rra_gates
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.err
#
# Module 04 step 3: cross-region interpretation gates.
#
# This one job spans all three regions by design -- the gates ARE the
# cross-region rule (H3K9me3/quiescent shared only if all three survive; LINE/L1
# multi-region if at least two). It cannot run per cell.
#
#   RRA_COHORT=AA RRA_RUN_IDS=id1,id2,id3 sbatch _h/step_3_gates.sh

V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

COHORT=${RRA_COHORT:?RRA_COHORT must be set}
RUN_IDS=${RRA_RUN_IDS:?RRA_RUN_IDS must be set}
# Analysis code is executed from the run's immutable snapshot under _m, not
# from the live _h/. V2_RUN_CODE is exported by the submit driver and points at
# {run_dir}/code/_h. Without this, editing _h/ while a run is queued silently
# changes what that run executes, and the snapshot and recorded git_commit
# attest to code that never ran. Falls back to live _h/ for a hand-run sbatch.
H_DIR="${V2_RUN_CODE:-${REPO_DIR}/04_repeat_repressive_architecture/_h}"

log_job_info
run_r "${H_DIR}/03_apply_gates.R" --cohort "$COHORT" --run-ids "$RUN_IDS"
