#!/bin/bash
# Submit step_4 (cis genotype extraction) for one arm x region.
#
#   cd 01_vmr_catalog/_m && mkdir -p logs
#   COHORT=AA REGION=dlpfc RUN_ID=vmrcat-AA-dlpfc-20260816 ../_h/submit_step_4.sh
#
# Why this wrapper exists: the cis window and the cohort's .pgen prefix are the
# same for every task in the array, but reading them costs two `conda run
# Rscript` startups (~50s) while the plink2 extraction itself takes ~2s. The
# first full DLPFC array therefore burned 141 CPU-hours resolving the same two
# values 19,144 times. Resolving them ONCE here and exporting them cuts the step
# to roughly 3s per task.
#
# It also sizes the array from vmr.bed rather than trusting a hand-typed N: an
# array shorter than the catalog silently skips VMRs, which would only surface
# much later as missing loci in 02_local_genetic_variance.

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT=AA|all_individuals}"
: "${REGION:?set REGION}"
: "${RUN_ID:?set RUN_ID}"

THROTTLE="${THROTTLE:-250}"

RUN_DIR="$REPO_DIR/01_vmr_catalog/_m/runs/$RUN_ID"
REGION_LIST="$RUN_DIR/vmr/vmr.bed"
require_file "$REGION_LIST"

N=$(wc -l < "$REGION_LIST")
if [ "$N" -lt 1 ]; then
    echo "ERROR: $REGION_LIST is empty" >&2
    exit 1
fi

config_value() {
    conda run --no-capture-output -p "$V2_ENV_R" Rscript -e \
        "source('$REPO_DIR/00_shared/load.R'); $1" | tail -1
}

V2_CIS_WINDOW_BP=$(config_value 'cat(load_config("thresholds")$cis$window_bp)')
: "${V2_CIS_WINDOW_BP:?could not read cis.window_bp from config/thresholds.yml}"

V2_PGEN_PREFIX=$(config_value "cat(cohort_def('$COHORT')\$pgen_prefix)")
: "${V2_PGEN_PREFIX:?could not resolve pgen_prefix for cohort $COHORT}"
case "$V2_PGEN_PREFIX" in
    /*) ;;
    *) V2_PGEN_PREFIX="$REPO_DIR/$V2_PGEN_PREFIX" ;;
esac
# Fail at submission, not on 9,572 tasks in sequence.
require_file "${V2_PGEN_PREFIX}.pgen"

export COHORT REGION RUN_ID V2_CIS_WINDOW_BP V2_PGEN_PREFIX V2_REPO_ROOT

echo "run       $RUN_ID"
echo "cohort    $COHORT / $REGION"
echo "VMRs      $N (throttle %${THROTTLE})"
echo "window    ${V2_CIS_WINDOW_BP} bp"
echo "pfile     ${V2_PGEN_PREFIX}"

exec sbatch --array=1-${N}%${THROTTLE} "$REPO_DIR/01_vmr_catalog/_h/step_4.sh"
