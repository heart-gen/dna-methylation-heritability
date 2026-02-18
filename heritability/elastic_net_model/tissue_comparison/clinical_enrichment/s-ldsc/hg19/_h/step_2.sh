#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=make_annot
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/make_annot.%j.log

# =============================================================================
# Step 2: Create LDSC Annotation Files
# =============================================================================
# This script creates annotation files for S-LDSC analysis using make_annot.py.
# Annotations are created for heritable and non-heritable VMRs in each brain
# region, for all 22 autosomes.
# =============================================================================

# Source configuration
PROJECT_BASE="/projects/b1213/users/elisa/dna-methylation-heritability"
SCRIPT_DIR="${PROJECT_BASE}/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

module purge
#module load anaconda3/2024.10-1
module load bedtools/2.31.1
module list

log_message "**** Loading conda environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

# Set temp directory for pybedtools (prevents /var/tmp issues on HPC)
# Options: Set USE_SHARED_TMP=1 to use shared project tmp, otherwise uses SCRATCH
USE_SHARED_TMP=1
SHARED_TMP="/projects/b1213/tmp"

if [[ "$USE_SHARED_TMP" -eq 1 ]] && [[ -d "$SHARED_TMP" ]]; then
    export TMPDIR="$SHARED_TMP"
elif [[ -n "$SCRATCH" ]] && [[ -d "$SCRATCH" ]]; then
    export TMPDIR="$SCRATCH"
else
    export TMPDIR="/tmp"
fi
mkdir -p "$TMPDIR"
log_message "Using TMPDIR: $TMPDIR"

# Validate resources
if ! validate_resources; then
    log_message "ERROR: Resource validation failed. Exiting."
    exit 1
fi

# LDSC annotation script
ANNOT_SCRIPT="${LDSC_DIR}/make_annot.py"

# Resource paths
BIM_DIR="${RESOURCE_DIR}/1000G_EUR_Phase3_plink"

# Process each brain region and heritability status
for REGION in "${BRAIN_REGIONS[@]}"; do
    log_message "Processing region: $REGION"

    for STATUS in "${HERITABILITY[@]}"; do
        log_message "  Processing status: $STATUS"

        # Input BED file
        BED_FILE="./vmr/${REGION}/${STATUS}.bed"

        # Output directory
        OUT_DIR="./custom_ldscores/${REGION}/${STATUS}"
        mkdir -p "$OUT_DIR"

        # Check if BED file exists
        if [[ ! -f "$BED_FILE" ]]; then
            log_message "    WARNING: BED file not found: $BED_FILE"
            continue
        fi

        # Filter and sort BED file (required for bedtools merge in make_annot.py)
        # Remove entries with empty chromosome (failed liftover regions)
        SORTED_BED_FILE="${OUT_DIR}/${REGION}_${STATUS}.sorted.bed"
        if [[ ! -f "$SORTED_BED_FILE" ]]; then
            log_message "    Filtering and sorting BED file..."
            awk -F'\t' '$1 != "" && $1 ~ /^chr[0-9]+$/' "$BED_FILE" | bedtools sort > "$SORTED_BED_FILE"
            log_message "    Filtered BED file: $(wc -l < "$SORTED_BED_FILE") valid entries"
        fi

        # Process each chromosome
        for CHR in {1..22}; do
            log_message "    Processing chromosome $CHR..."

            BIM_FILE="${BIM_DIR}/1000G.EUR.QC.${CHR}.bim"
            ANNOT_FILE="${OUT_DIR}/${REGION}_${STATUS}.${CHR}.annot.gz"

            # Check if BIM file exists
            if [[ ! -f "$BIM_FILE" ]]; then
                log_message "      WARNING: BIM file not found: $BIM_FILE"
                continue
            fi

            # Skip if annotation already exists
            if [[ -f "$ANNOT_FILE" ]]; then
                log_message "      Annotation already exists, skipping..."
                continue
            fi

            python "$ANNOT_SCRIPT" \
                --bed-file "$SORTED_BED_FILE" \
                --bimfile "$BIM_FILE" \
                --annot-file "$ANNOT_FILE" \
                --windowsize 500000
        done
    done
done

# -----------------------------------------------------------------------------
# Verify outputs
# -----------------------------------------------------------------------------
log_message "Verifying annotation files..."
echo ""
echo "=== Annotation Files ==="
for REGION in "${BRAIN_REGIONS[@]}"; do
    for STATUS in "${HERITABILITY[@]}"; do
        OUT_DIR="./custom_ldscores/${REGION}/${STATUS}"
        count=$(ls -1 "${OUT_DIR}"/*.annot.gz 2>/dev/null | wc -l)
        if [[ $count -eq 22 ]]; then
            echo "[OK] ${REGION}/${STATUS}: $count annotation files"
        else
            echo "[INCOMPLETE] ${REGION}/${STATUS}: $count/22 annotation files"
        fi
    done
done

log_message "**** Job ends ****"
