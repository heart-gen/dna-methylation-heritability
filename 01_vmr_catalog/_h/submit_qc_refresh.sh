#!/bin/bash
#
# 01_vmr_catalog: refresh QC tables for every accepted cohort x region.
#
# Each cell is independent, so all six are submitted at once. They mint their
# own run IDs and seal themselves; there is nothing to chain.
#
# Usage, from the module's _m directory:
#   cd 01_vmr_catalog/_m && mkdir -p logs
#   ../_h/submit_qc_refresh.sh                    # all six cells
#   SOURCE_DATE=20260816 ../_h/submit_qc_refresh.sh
#
# Environment:
#   SOURCE_DATE=YYYYMMDD  date stamp of the accepted catalog runs (default 20260816)
#   DRY_RUN=1             print the plan without submitting

set -euo pipefail

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

SOURCE_DATE="${SOURCE_DATE:-20260816}"
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

for COHORT in AA all_individuals; do
    for REGION in caudate dlpfc hippocampus; do
        SOURCE_RUN_ID="vmrcat-${COHORT}-${REGION}-${SOURCE_DATE}"
        if [ ! -d "$REPO_DIR/01_vmr_catalog/_m/runs/$SOURCE_RUN_ID" ]; then
            echo "ERROR: no such source run: $SOURCE_RUN_ID" >&2
            exit 1
        fi
        JOB=$(submit --export="ALL,COHORT=$COHORT,REGION=$REGION,SOURCE_RUN_ID=$SOURCE_RUN_ID" \
              "$HERE/step_4b_qc_refresh.sh")
        log_message "${COHORT}/${REGION} from ${SOURCE_RUN_ID}: $JOB"
    done
done

cat <<EOF

Submitted six QC refresh runs.

Each mints its own run ID (vmrcatqc-{cohort}-{region}-{date}) and seals itself.
Collect the IDs from the job logs, then point the Figure 1 builder at them:
  10_integrated_manuscript_outputs/_h/01_figure1_catalog.R  (QC_RUN)
EOF
