#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=manuscript_figures
#SBATCH --output=logs/manuscript_figures.%A.log
#
# 10_integrated_manuscript_outputs step 1: assemble Figures 1 and 2.
#
# Consumes only accepted immutable upstream runs (AGENTS.md 7.9). Figure 1
# needs the Module 01 catalog runs AND the QC refresh runs that carry array
# coverage and genomic context; Figure 2 needs the Module 02 score runs.
#
# The builders resolve those upstream run IDs internally, so the only argument
# is the output run ID.
#
# Usage, from the module's _m directory:
#   cd 10_integrated_manuscript_outputs/_m && mkdir -p logs
#   RUN_ID=fig-all-20260826-a sbatch ../_h/step_1_figures.sh

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${RUN_ID:?set RUN_ID}"

HERE="$REPO_DIR/10_integrated_manuscript_outputs/_h"
mkdir -p "$REPO_DIR/10_integrated_manuscript_outputs/_m/runs/$RUN_ID"

log_job_info
log_message "**** Building manuscript figures into ${RUN_ID} ****"

# AA is the primary arm; all_individuals renders from the same builders as the
# sensitivity supplement. Figure 1 gets a second pass with EPIC as the stricter
# array comparator.
for COHORT in AA all_individuals; do
    log_message "Figure 1 (450K) -- ${COHORT}"
    run_r "$HERE/01_figure1_catalog.R" --cohort "$COHORT" --run-id "$RUN_ID"

    log_message "Figure 1 (EPIC supplement) -- ${COHORT}"
    run_r "$HERE/01_figure1_catalog.R" --cohort "$COHORT" --run-id "$RUN_ID" \
        --platform EPIC

    log_message "Figure 2 + audit supplement -- ${COHORT}"
    run_r "$HERE/02_figure2_local_control.R" --cohort "$COHORT" --run-id "$RUN_ID"
done

# Writes the provenance manifest, verifies every figure has source data, and
# seals the run read-only. Must be last.
log_message "Sealing ${RUN_ID}"
run_r "$HERE/03_close_figure_run.R" --run-id "$RUN_ID"

log_message "**** Job ends ****"
