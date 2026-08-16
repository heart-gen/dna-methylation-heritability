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
#
# load.R is addressed absolutely: this step is submitted from _m/, which is the
# module convention, so a relative "00_shared/load.R" does not exist from the
# job's working directory. Stderr is NOT discarded -- with `set -e`, a failing
# command substitution kills the script at the assignment, before the :? guard
# below can report anything, so suppressing it produced a completely empty log.
config_value() {
    conda run --no-capture-output -p "$V2_ENV_R" Rscript -e \
        "source('$REPO_DIR/00_shared/load.R'); $1" | tail -1
}

WINDOW=$(config_value 'cat(load_config("thresholds")$cis$window_bp)')
: "${WINDOW:?could not read cis.window_bp from config/thresholds.yml}"

PFILE=$(config_value "cat(cohort_def('$COHORT')\$pgen_prefix)")
: "${PFILE:?could not resolve pgen_prefix for cohort $COHORT}"
# cohort_def() already returns an absolute path; only prefix REPO_DIR if some
# future config makes it repo-relative, rather than doubling it unconditionally.
case "$PFILE" in
    /*) ;;
    *) PFILE="$REPO_DIR/$PFILE" ;;
esac
require_file "${PFILE}.pgen"

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

# plink2 defaults to the node's core count, not the allocation, so it needs the
# ceiling passed explicitly -- same defect as data.table (00_shared/slurm.sh).
OUT_PREFIX="$CHR_DIR/${COHORT}.${START}_${END}"
STATUS="extracted"

# plink2 defaults to the node's core count, not the allocation, so it needs the
# ceiling passed explicitly -- same defect as data.table (00_shared/slurm.sh).
PLINK_RC=0
"$PLINK2" --pfile "$PFILE" \
          --chr "${CHR#chr}" \
          --from-bp "$START_POS" \
          --to-bp "$END_POS" \
          --keep "$RUN_DIR/vmr/donors_plink.txt" \
          --make-bed \
          --no-parents \
          --no-sex \
          --no-pheno \
          --out "$OUT_PREFIX" || PLINK_RC=$?

# A VMR whose cis window contains no genotyped variants is an EXPLAINED
# EXCLUSION, not a task failure. It happens for real: the first genotyped
# variant on chr1 is at 833,068, so the four DLPFC VMRs below that position have
# an empty window however wide it is. Such a VMR simply has no local genetic
# variance to estimate, and 02_local_genetic_variance must skip it.
#
# Failing the task instead would put a permanent non-zero `failed` count into
# task_reconciliation.tsv, which AGENTS.md 6 treats as blocking -- an unexplained
# failure and a locus with no cis SNPs would become indistinguishable.
if [ "$PLINK_RC" -ne 0 ]; then
    if grep -q "No variants remaining after main filters" "${OUT_PREFIX}.log" 2>/dev/null; then
        STATUS="no_cis_variants"
        log_message "no genotyped variants in ${CHR}:${START_POS}-${END_POS}; recording exclusion"
        rm -f "${OUT_PREFIX}.bed" "${OUT_PREFIX}.bim" "${OUT_PREFIX}.fam"
    else
        echo "ERROR: plink2 failed (exit $PLINK_RC) for ${CHR}:${START}-${END}" >&2
        exit "$PLINK_RC"
    fi
else
    require_file "${OUT_PREFIX}.bed"
    N_VAR=$(wc -l < "${OUT_PREFIX}.bim")
    log_message "extracted ${N_VAR} variants"
fi

# Record what was actually extracted, so the combine step can reconcile.
# One file per task, not one shared append -- concurrent array tasks appending to
# a single file interleave and lose lines.
mkdir -p "$OUTPUT/extraction_log"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$TASK" "$CHR" "$START" "$END" "$START_POS" "$END_POS" "$CLAMPED" \
    "$STATUS" "${N_VAR:-0}" \
    > "$OUTPUT/extraction_log/task_${TASK}.tsv"

log_message "**** Job ends ****"
