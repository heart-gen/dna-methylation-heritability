#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name=make_annot
#SBATCH --output=logs/02.output_%j.log
#SBATCH --error=logs/02.error_%j.log

# =============================================================================
# Step 2: Create LDSC Annotation Files
# =============================================================================
# This script creates annotation files for S-LDSC analysis using make_annot.py.
# Annotations are created for heritable and non-heritable VMRs in each brain
# region, for all 22 autosomes.
# =============================================================================

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

# Load bedtools module
module load bedtools/2.30.0 2>/dev/null || true

# Validate resources
if ! validate_resources; then
    log_message "ERROR: Resource validation failed. Exiting."
    exit 1
fi

# LDSC annotation script
ANNOT_SCRIPT="${LDSC_DIR}/make_annot.py"

# Resource paths
BIM_DIR="${RESOURCE_DIR}/1000G_Phase3_plinkfiles"

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
                --bed-file "$BED_FILE" \
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
