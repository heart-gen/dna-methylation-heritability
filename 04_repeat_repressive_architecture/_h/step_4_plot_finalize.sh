#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=rra_plot_final
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.err
#
# Module 04 step 4: figures, then seal every cell.
#
# Like the gates, this spans all three regions: the figures are drawn from the
# pooled cross-region tables, and no cell may be sealed until the gates that
# decide what it is allowed to claim have been applied to all three.
#
#   RRA_COHORT=AA RRA_RUN_IDS=id1,id2,id3 sbatch _h/step_4_plot_finalize.sh

V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

COHORT=${RRA_COHORT:?RRA_COHORT must be set}
RUN_IDS=${RRA_RUN_IDS:?RRA_RUN_IDS must be set}
H_DIR="${REPO_DIR}/04_repeat_repressive_architecture/_h"

log_job_info
run_r "${H_DIR}/04_plot.R" --cohort "$COHORT" --run-ids "$RUN_IDS"
run_r "${H_DIR}/05_finalize_run.R" --cohort "$COHORT" --run-ids "$RUN_IDS"
