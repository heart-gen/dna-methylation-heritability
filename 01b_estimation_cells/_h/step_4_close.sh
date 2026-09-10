#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=estcell_close
#SBATCH --output=logs/estcell_close.%j.log
#
# Takes no positional arguments; reads RUN_ID from the environment, so any
# single stage can be resubmitted by hand against an existing run without
# editing the workflow script.
#
#   RUN_ID=<id> sbatch _h/step_4_close.sh

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${RUN_ID:?set RUN_ID}"

log_job_info
run_r "$REPO_DIR/01b_estimation_cells/_h/04_close_run.R" --run-id "$RUN_ID"
log_message "**** Job ends ****"
