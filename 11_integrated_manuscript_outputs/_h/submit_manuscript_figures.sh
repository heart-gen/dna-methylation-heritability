#!/bin/bash
#
# 11_integrated_manuscript_outputs: mint a figure run ID and submit the build.
#
# Usage, from the module's _m directory:
#   cd 11_integrated_manuscript_outputs/_m && mkdir -p logs
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

HERE="$REPO_DIR/11_integrated_manuscript_outputs/_h"
RUNS="$REPO_DIR/11_integrated_manuscript_outputs/_m/runs"
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

# Fail before queueing if any upstream module the builders consume lacks an
# accepted run for a cell they need.
#
# This used to require_file() six hard-coded lgv-*-20260823 paths. Module 02
# retired those on 2026-09-17 in favour of the rescored runs, and because the
# check named files rather than asking the acceptance gate, it kept passing
# against a superseded run -- which is how the last build of Figure 2 shipped
# uncitable. The gate is the record of decision (AGENTS.md 6); ask it.
log_message "Preflight: resolving accepted upstream runs"
conda run --no-capture-output -p "$V2_ENV_R" Rscript - <<'RS'
suppressPackageStartupMessages(source(file.path(Sys.getenv("V2_REPO_ROOT"), "00_shared", "load.R")))
regions <- load_config("cohorts")$regions

## Module -> the cohorts whose cells the build needs. Figures 1-2 render for
## both arms; Figures 3-5 and the supplements are AA-only, because those
## modules have accepted runs for AA alone.
need <- list(
    "01_vmr_catalog"                     = c("AA", "all_individuals"),
    "02_local_genetic_variance"          = c("AA", "all_individuals"),
    "04_repeat_repressive_architecture"  = "AA",
    "05_cpg_meqtl_burden"                = "AA",
    "06_partitioned_heritability"        = "AA",
    "07_transcription_splicing_coupling" = "AA",
    "09_schizophrenia_risk_application"  = "AA",
    "09b_aging_application"              = "AA",
    "10_environmental_exploratory"       = "AA")

bad <- character(0)
for (m in names(need)) for (co in need[[m]]) for (re in regions) {
    r <- tryCatch(require_accepted_upstream(m, co, re)$run_id,
                  error = function(e) NA_character_)
    if (is.na(r)) bad <- c(bad, sprintf("%s [%s x %s]", m, co, re))
    else message(sprintf("  %-36s %-16s %-12s %s", m, co, re, r))
}

## Module 08 spans all three regions under the literal region `crossregion`.
r8 <- tryCatch(require_accepted_upstream("08_region_donor_generalization", "AA",
                                         "crossregion")$run_id,
               error = function(e) NA_character_)
## Braced: at top level R closes `if (cond) expr` before it sees `else`.
if (is.na(r8)) {
    bad <- c(bad, "08_region_donor_generalization [AA x crossregion]")
} else {
    message(sprintf("  %-36s %-16s %-12s %s", "08_region_donor_generalization",
                    "AA", "crossregion", r8))
}

## Figure 5 reads the module-level stage 17/18 tables. They carry -UNACCEPTED
## until both stages are re-run without --allow-unlocked.
sfx <- file.path(repo_root(), "09_schizophrenia_risk_application", "_m", "combined")
for (f in c("scz-negative-control-summary-AA.tsv",
            "scz-locus-architecture-axis-link-AA.tsv")) {
    if (!file.exists(file.path(sfx, f))) {
        bad <- c(bad, paste0("09 stage 17/18: ", f, " (run _h/step_9_negative_controls.sh)"))
    }
}

if (length(bad) > 0) {
    stop("No accepted run for:\n  ", paste(bad, collapse = "\n  "),
         "\nRecord the run ID and its passing gate in the module README ",
         "before building manuscript figures on it (AGENTS.md 6).")
}
message("Preflight OK")
RS

log_message "RUN_ID=${RUN_ID}"
if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "[dry-run] sbatch --export=ALL,RUN_ID=$RUN_ID $HERE/step_1_figures.sh" >&2
    exit 0
fi

# Snapshot the code and configuration into the run, and execute the snapshot
# rather than the live tree -- the convention Modules 04, 09 and 10 already
# follow (09/_h/submit_schizophrenia.sh:39-47).
#
# This module lacked it, and that is what made fig-all-20260920 unacceptable:
# its manifest recorded git_dirty = true against a commit that did not contain
# the builders, so the only record of the code that produced the run was an
# uncommitted working tree. A figure run is the one place where that matters
# most, because every manuscript number is traced back through it.
#
# Snapshotting also closes the narrower race the other modules hit: editing _h/
# or config/ while this job is queued would otherwise change what the run
# executes while its manifest still attests to the original checksums.
RUN_DIR="$RUNS/$RUN_ID"
mkdir -p "$RUN_DIR/code"
cp -a "$HERE" "$RUN_DIR/code/_h"
mkdir -p "$RUN_DIR/code/config"
cp -a "$REPO_DIR/config/." "$RUN_DIR/code/config/"
RUN_CODE="$RUN_DIR/code/_h"
log_message "snapshotted _h/ and config/ into ${RUN_DIR}/code"

JOB=$(sbatch --parsable --export="ALL,RUN_ID=$RUN_ID,V2_RUN_CODE=$RUN_CODE" \
    "$RUN_CODE/step_1_figures.sh")
log_message "step_1 (figures + seal): $JOB"

cat <<EOF

Submitted manuscript figures as run ${RUN_ID}

  run dir : $RUNS/$RUN_ID
  job     : $JOB

Figures land in figures/, per-panel provenance in source_data/. The run seals
itself; every panel records its source run ID, table, script, and filter.
EOF
