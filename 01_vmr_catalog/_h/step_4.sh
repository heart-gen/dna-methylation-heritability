#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=vmr_extract_snp
#SBATCH --output=logs/vmr_extract_snp.%A_%a.log
#
# 01_vmr_catalog step 4: extract the cis genotype window for each VMR.
#
# Replaces the legacy vmr-analysis/*/_h/step_5.sh. Repairs:
#
#   V9  Legacy ABORTED when a VMR sat within one window of a chromosome start
#       (START_POS <= 0), silently dropping ~1% of VMRs -- 107/102/100 in
#       BA_only and 336/284/244 in all_individuals. v2 CLAMPS to 1 instead.
#   V10 Legacy looked up chr1's length once and used it as the upper bound for
#       every chromosome. chr1 is the longest, so the bound never bit on smaller
#       chromosomes and the check was effectively dead. v2 looks up each
#       chromosome via chrom_size() and clamps to it.
#   V11 The window came from a literal that differed between arms (500 kb in
#       BA_only, 1 Mb in all_individuals), making the two cohorts incomparable.
#       v2 reads it from config/thresholds.yml, so both arms share one value.
#   V12 plink2 from the shared opt tree, never `module load`.
#
# The array size must match the VMR count for this run:
#   N=$(wc -l < .../vmr/vmr.bed)
#   sbatch --array=1-${N}%250 ../_h/step_4.sh

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${COHORT:?set COHORT=AA|all_individuals}"
: "${REGION:?set REGION}"
: "${RUN_ID:?set RUN_ID}"

RUN_DIR="$REPO_DIR/01_vmr_catalog/_m/runs/$RUN_ID"
REGION_LIST="$RUN_DIR/vmr/vmr.bed"
OUTPUT="$RUN_DIR/plink_format"

require_file "$REGION_LIST"
require_exec "$PLINK2"

# Window and cohort genotype prefix come from config, not from literals.
WINDOW=$(conda run --no-capture-output -p "$V2_ENV_R" Rscript -e \
    'source("00_shared/load.R"); cat(load_config("thresholds")$cis$window_bp)' \
    2>/dev/null | tail -1)
: "${WINDOW:?could not read cis.window_bp from config/thresholds.yml}"

PFILE=$(conda run --no-capture-output -p "$V2_ENV_R" Rscript -e \
    "source(\"00_shared/load.R\"); cat(cohort_def(\"$COHORT\")\$pgen_prefix)" \
    2>/dev/null | tail -1)
: "${PFILE:?could not resolve pgen_prefix for cohort $COHORT}"

TASK="${SLURM_ARRAY_TASK_ID:?this step must run as a SLURM array}"
LINE=$(sed -n "${TASK}p" "$REGION_LIST")
if [ -z "$LINE" ]; then
    echo "ERROR: no VMR on line $TASK of $REGION_LIST" >&2
    exit 1
fi

CHR=$(echo "$LINE" | awk '{print $1}')
START=$(echo "$LINE" | awk '{print $2}')
END=$(echo "$LINE" | awk '{print $3}')

log_job_info
log_message "VMR ${CHR}:${START}-${END} (task ${TASK}, window ${WINDOW} bp)"

# V10: this chromosome's length, not chr1's.
CHR_SIZE=$(chrom_size "$CHR")

# V9: clamp rather than abort. A VMR near a telomere gets a truncated window,
# which is a smaller cis region -- not a reason to drop the VMR entirely.
START_POS=$((START - WINDOW))
END_POS=$((END + WINDOW))
CLAMPED="none"
if (( START_POS < 1 )); then
    START_POS=1
    CLAMPED="start"
fi
if (( END_POS > CHR_SIZE )); then
    END_POS=$CHR_SIZE
    CLAMPED="${CLAMPED}+end"
fi
if [ "$CLAMPED" != "none" ]; then
    log_message "window clamped (${CLAMPED}) to ${START_POS}-${END_POS} of ${CHR} (${CHR_SIZE} bp)"
fi

CHR_DIR="$OUTPUT/${CHR}"
mkdir -p "$CHR_DIR"

"$PLINK2" --pfile "$REPO_DIR/$PFILE" \
          --chr "${CHR#chr}" \
          --from-bp "$START_POS" \
          --to-bp "$END_POS" \
          --keep "$RUN_DIR/vmr/donors_plink.txt" \
          --make-bed \
          --no-parents \
          --no-sex \
          --no-pheno \
          --out "$CHR_DIR/${COHORT}.${START}_${END}"

# Record what was actually extracted, so the combine step can reconcile.
# One file per task, not one shared append -- concurrent array tasks appending to
# a single file interleave and lose lines.
mkdir -p "$OUTPUT/extraction_log"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$TASK" "$CHR" "$START" "$END" "$START_POS" "$END_POS" "$CLAMPED" \
    > "$OUTPUT/extraction_log/task_${TASK}.tsv"

log_message "**** Job ends ****"
