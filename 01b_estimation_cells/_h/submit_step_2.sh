#!/bin/bash
# Submit step 2 (per-cell cis genotype extraction) for one estimation cell.
#
#   cd 01b_estimation_cells/_m && mkdir -p logs
#   RUN_ID=estcell-all_individuals.EA-dlpfc-20260910 ../_h/submit_step_2.sh
#
# Mirrors 01_vmr_catalog/_h/submit_step_4.sh: resolve the values that are
# identical for every array task ONCE and export them, and size the array from
# vmr.bed rather than a hand-typed N. An array shorter than the catalog silently
# skips VMRs, which would only surface much later as missing loci in Module 02.
#
# Unlike step_4.sh's wrapper these values come from the run manifest, not from
# config: the run already recorded which pgen and which window it was opened
# with, and a run must not change its inputs after the fact.

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${RUN_ID:?set RUN_ID}"

THROTTLE="${THROTTLE:-250}"

RUN_DIR="$REPO_DIR/01b_estimation_cells/_m/runs/$RUN_ID"
REGION_LIST="$RUN_DIR/vmr/vmr.bed"
require_file "$REGION_LIST"
require_file "$RUN_DIR/manifest.tsv"
require_file "$RUN_DIR/vmr/donors_plink.txt"

N=$(wc -l < "$REGION_LIST")
if [ "$N" -lt 1 ]; then
    echo "ERROR: $REGION_LIST is empty" >&2
    exit 1
fi

manifest_value() {
    awk -F'\t' -v f="$1" '$1 == f { print $2; found = 1 } END { if (!found) exit 1 }' \
        "$RUN_DIR/manifest.tsv"
}

V2_CIS_WINDOW_BP=$(manifest_value cis_window_bp)
V2_PGEN_PREFIX=$(manifest_value pgen_prefix)
V2_ESTIMATION_GROUP=$(manifest_value estimation_group)
: "${V2_CIS_WINDOW_BP:?manifest carries no cis_window_bp}"
: "${V2_PGEN_PREFIX:?manifest carries no pgen_prefix}"
: "${V2_ESTIMATION_GROUP:?manifest carries no estimation_group}"

case "$V2_PGEN_PREFIX" in
    /*) ;;
    *) V2_PGEN_PREFIX="$REPO_DIR/$V2_PGEN_PREFIX" ;;
esac
# Fail at submission, not on N tasks in sequence.
require_file "${V2_PGEN_PREFIX}.pgen"

export RUN_ID V2_CIS_WINDOW_BP V2_PGEN_PREFIX V2_ESTIMATION_GROUP V2_REPO_ROOT

echo "run       $RUN_ID"
echo "cell      $(manifest_value cohort) / $(manifest_value region)"
echo "group     $V2_ESTIMATION_GROUP ($(manifest_value n_donors_cell) of $(manifest_value n_donors_pooled) donors)"
echo "VMRs      $N (throttle %${THROTTLE})"
echo "window    ${V2_CIS_WINDOW_BP} bp"
echo "pfile     ${V2_PGEN_PREFIX}"

exec sbatch --array=1-${N}%${THROTTLE} \
    "$REPO_DIR/01b_estimation_cells/_h/step_2_extract_cis.sh"
