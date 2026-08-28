#!/bin/bash
#
# 10_integrated_manuscript_outputs: mint a figure run ID and submit the build.
#
# Usage, from the module's _m directory:
#   cd 10_integrated_manuscript_outputs/_m && mkdir -p logs
#   ../_h/submit_manuscript_figures.sh
#   RUN_ID=fig-all-20260826-b ../_h/submit_manuscript_figures.sh
#
# Environment:
#   RUN_ID=...  use this run ID instead of minting one
#   DRY_RUN=1   print the plan without submitting

set -euo pipefail

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

HERE="$REPO_DIR/10_integrated_manuscript_outputs/_h"
RUNS="$REPO_DIR/10_integrated_manuscript_outputs/_m/runs"
mkdir -p logs

# Run directories are immutable, so never reuse one. Suffix until free, the
# same rule make_run_id() applies in 00_shared/runid.R.
if [ -z "${RUN_ID:-}" ]; then
    BASE="fig-all-$(date +%Y%m%d)"
    RUN_ID="$BASE"
    for SUFFIX in a b c d e f g h i j; do
        [ ! -d "$RUNS/$RUN_ID" ] && break
        RUN_ID="${BASE}-${SUFFIX}"
    done
fi
if [ -d "$RUNS/$RUN_ID" ]; then
    echo "ERROR: run directory already exists and runs are immutable: $RUNS/$RUN_ID" >&2
    exit 1
fi

# Fail before queueing if an upstream run the builders need is absent.
for COHORT in AA all_individuals; do
    for REGION in caudate dlpfc hippocampus; do
        require_file "$REPO_DIR/02_local_genetic_variance/_m/runs/lgv-${COHORT}-${REGION}-20260823/results/combined/local-genetic-control-${COHORT}-${REGION}-vmrs.tsv"
    done
done

log_message "RUN_ID=${RUN_ID}"
if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "[dry-run] sbatch --export=ALL,RUN_ID=$RUN_ID $HERE/step_1_figures.sh" >&2
    exit 0
fi

JOB=$(sbatch --parsable --export="ALL,RUN_ID=$RUN_ID" "$HERE/step_1_figures.sh")
log_message "step_1 (figures + seal): $JOB"

cat <<EOF

Submitted manuscript figures as run ${RUN_ID}

  run dir : $RUNS/$RUN_ID
  job     : $JOB

Figures land in figures/, per-panel provenance in source_data/. The run seals
itself; every panel records its source run ID, table, script, and filter.
EOF
