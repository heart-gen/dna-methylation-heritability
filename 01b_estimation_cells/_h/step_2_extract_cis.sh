#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=estcell_extract
#SBATCH --output=logs/estcell_extract.%A_%a.log
#
# 01b step 2: extract each VMR's cis genotype window for THIS cell's donors.
#
# Same windows, same pooled pgen, same repairs as 01_vmr_catalog/_h/step_4.sh,
# which this is adapted from and must stay behaviourally identical to:
#
#   V9  clamp START_POS to 1 instead of aborting near a chromosome start
#   V10 look up THIS chromosome's length, not chr1's
#   V11 window from config/thresholds.yml, never a literal
#   V12 plink2 from the shared opt tree, never `module load`
#
# The only difference from step_4.sh is --keep: this cell's donor list rather
# than the pooled one. That is deliberately the ONLY difference. Downstream MAF
# and missingness QC then run on the group-restricted matrix inside
# 00_shared/locus_io.R, which is how the >= maf_min filter becomes within-group.
#
# This is the v1 design (vmr-analysis/all_individuals/{region}/_h/step_5.sh:95-122
# ran plink2 twice per window with --keep samples-AA.txt / samples-EA.txt),
# rebuilt on the v2 run contract.
#
# Submit via the wrapper, which sizes the array and resolves config once:
#   cd 01b_estimation_cells/_m && mkdir -p logs
#   RUN_ID=<id> ../_h/submit_step_2.sh

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${RUN_ID:?set RUN_ID}"

RUN_DIR="$REPO_DIR/01b_estimation_cells/_m/runs/$RUN_ID"
REGION_LIST="$RUN_DIR/vmr/vmr.bed"
KEEP="$RUN_DIR/vmr/donors_plink.txt"
OUTPUT="$RUN_DIR/plink_format"

require_file "$REGION_LIST"
require_file "$KEEP"
require_exec "$PLINK2"

manifest_value() {
    awk -F'\t' -v f="$1" '$1 == f { print $2; found = 1 } END { if (!found) exit 1 }' \
        "$RUN_DIR/manifest.tsv"
}

# Resolved once by submit_step_2.sh and inherited; the lookup is the fallback
# for a task run standalone. See step_4.sh for why this matters at array scale.
WINDOW="${V2_CIS_WINDOW_BP:-$(manifest_value cis_window_bp)}"
: "${WINDOW:?could not resolve cis_window_bp}"

PFILE="${V2_PGEN_PREFIX:-$(manifest_value pgen_prefix)}"
: "${PFILE:?could not resolve pgen_prefix}"
require_file "${PFILE}.pgen"

# The BED filename token is the ESTIMATION GROUP, which is what
# 00_shared/locus_io.R builds its path from. AA / EA never collide with the
# pooled arm's own all_individuals files, so a cell's extraction is unambiguous.
GROUP="${V2_ESTIMATION_GROUP:-$(manifest_value estimation_group)}"
: "${GROUP:?could not resolve estimation_group}"

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
log_message "VMR ${CHR}:${START}-${END} (task ${TASK}, group ${GROUP}, window ${WINDOW} bp)"

# V10: this chromosome's length, not chr1's.
CHR_SIZE=$(chrom_size "$CHR")

# V9: clamp rather than abort.
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

CHR_DIR="$OUTPUT/chr_${CHR#chr}"
mkdir -p "$CHR_DIR"

OUT_PREFIX="$CHR_DIR/TOPMed_LIBD-${GROUP}.${START}_${END}"
STATUS="extracted"

PLINK_RC=0
"$PLINK2" --pfile "$PFILE" \
          --chr "${CHR#chr}" \
          --from-bp "$START_POS" \
          --to-bp "$END_POS" \
          --keep "$KEEP" \
          --make-bed \
          --no-parents \
          --no-sex \
          --no-pheno \
          --out "$OUT_PREFIX" || PLINK_RC=$?

# A window with no genotyped variants is an EXPLAINED EXCLUSION, not a failure.
# Note this can differ BETWEEN cells on the same catalog: a window can carry
# variants in one donor group and none in the other. That asymmetry is real and
# must be recorded per cell, which is why the marker is written here rather than
# inherited from the pooled run.
if [ "$PLINK_RC" -ne 0 ]; then
    if grep -q "No variants remaining after main filters" "${OUT_PREFIX}.log" 2>/dev/null; then
        STATUS="no_cis_variants"
        log_message "no genotyped variants in ${CHR}:${START_POS}-${END_POS}; recording exclusion"
        rm -f "${OUT_PREFIX}.bed" "${OUT_PREFIX}.bim" "${OUT_PREFIX}.fam"
        : > "${OUT_PREFIX}.no-snps"
    else
        echo "ERROR: plink2 failed (exit $PLINK_RC) for ${CHR}:${START}-${END}" >&2
        exit "$PLINK_RC"
    fi
else
    require_file "${OUT_PREFIX}.bed"
    N_VAR=$(wc -l < "${OUT_PREFIX}.bim")
    log_message "extracted ${N_VAR} variants"
fi

# One file per task, never a shared append: concurrent array tasks interleave.
mkdir -p "$OUTPUT/extraction_log"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$TASK" "$CHR" "$START" "$END" "$START_POS" "$END_POS" "$CLAMPED" \
    "$STATUS" "${N_VAR:-0}" \
    > "$OUTPUT/extraction_log/task_${TASK}.tsv"

log_message "**** Job ends ****"
