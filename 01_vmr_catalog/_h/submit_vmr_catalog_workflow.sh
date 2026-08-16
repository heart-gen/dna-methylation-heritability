#!/bin/bash
#
# 01_vmr_catalog: submit the full catalog for one cohort x region, with SLURM
# dependencies chaining the steps.
#
# Usage, from the module's _m/ directory:
#   cd 01_vmr_catalog/_m && mkdir -p logs
#   ../_h/submit_vmr_catalog_workflow.sh AA caudate
#
# Environment:
#   ALLOW_UNLOCKED=1   smoke run; bypasses the pi_locked gate. Output must not
#                      be used as production or cited downstream.
#   DRY_RUN=1          print the plan without submitting.
#   WITH_SEX=1         also prepare X/Y into the run's excluded/ directory.
#   WITH_GENOTYPES=1   also run step_4 genotype extraction after step_3.

set -euo pipefail

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

COHORT="${1:?usage: submit_vmr_catalog_workflow.sh <AA|all_individuals> <region>}"
REGION="${2:?usage: submit_vmr_catalog_workflow.sh <cohort> <caudate|dlpfc|hippocampus>}"

HERE="$REPO_DIR/01_vmr_catalog/_h"
mkdir -p logs

submit() {
    if [ "${DRY_RUN:-0}" = "1" ]; then
        echo "[dry-run] sbatch $*" >&2
        echo "000000"
    else
        sbatch --parsable "$@"
    fi
}

log_message "Creating run for ${COHORT}/${REGION}"
RUN_ID=$(cd "$REPO_DIR" && conda run --no-capture-output -p "$V2_ENV_R" Rscript \
    "$HERE/00_new_run.R" --cohort "$COHORT" --region "$REGION" \
    ${ALLOW_UNLOCKED:+--allow-unlocked} | tail -1)
: "${RUN_ID:?run creation failed}"
log_message "RUN_ID=${RUN_ID}"

export COHORT REGION RUN_ID
[ "${ALLOW_UNLOCKED:-0}" = "1" ] && export ALLOW_UNLOCKED=1

EXPORTS="ALL,COHORT=$COHORT,REGION=$REGION,RUN_ID=$RUN_ID"
[ "${ALLOW_UNLOCKED:-0}" = "1" ] && EXPORTS="$EXPORTS,ALLOW_UNLOCKED=1"

J1=$(submit --export="$EXPORTS" "$HERE/step_1.sh")
log_message "step_1 (prepare, chr1-22): $J1"

DEPS="afterok:$J1"
if [ "${WITH_SEX:-0}" = "1" ]; then
    J1X=$(submit --export="$EXPORTS" "$HERE/step_1x.sh")
    log_message "step_1x (prepare X/Y -> excluded/): $J1X"
    DEPS="$DEPS,afterok:$J1X"
fi

J2=$(submit --dependency="$DEPS" --export="$EXPORTS" "$HERE/step_2.sh")
log_message "step_2 (PCs + residual variance): $J2"

J3=$(submit --dependency="afterok:$J2" --export="$EXPORTS" "$HERE/step_3.sh")
log_message "step_3 (call VMRs, plots, turnover): $J3"

if [ "${WITH_GENOTYPES:-0}" = "1" ]; then
    # step_4's array size depends on the VMR count, which does not exist until
    # step_3 finishes. Submit a small job that sizes and submits the array.
    log_message "step_4 (genotype extraction) is gated on the VMR count."
    log_message "After $J3 completes, run:"
    log_message "  N=\$(wc -l < $REPO_DIR/01_vmr_catalog/_m/runs/$RUN_ID/vmr/vmr.bed)"
    log_message "  COHORT=$COHORT REGION=$REGION RUN_ID=$RUN_ID \\"
    log_message "    sbatch --array=1-\${N}%250 $HERE/step_4.sh"
    log_message "Then seal the run (must be last -- it makes the run read-only):"
    log_message "  COHORT=$COHORT REGION=$REGION RUN_ID=$RUN_ID \\"
    log_message "    sbatch --dependency=afterok:<step4_jobid> $HERE/step_5.sh"
else
    J5=$(submit --dependency="afterok:$J3" --export="$EXPORTS" "$HERE/step_5.sh")
    log_message "step_5 (seal run): $J5"
fi

cat <<EOF

Submitted ${COHORT}/${REGION} as run ${RUN_ID}

  run dir : $REPO_DIR/01_vmr_catalog/_m/runs/$RUN_ID
  jobs    : prepare=$J1 analyze=$J2 summarize=$J3

Before treating this run as production, confirm in the module README:
  - observed donor counts match config/cohorts.yml design counts
  - task reconciliation reports zero unaccounted or failed tasks
  - turnover against the legacy catalog is reported and reviewed
EOF
