#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=GPU-shared
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --array=1-22
#SBATCH --gpus=1
#SBATCH --job-name=compute_ldscores
#SBATCH --output=logs/03.output_%a.log
#SBATCH --error=logs/03.error_%a.log

# =============================================================================
# Step 3: Compute LD Scores
# =============================================================================
# This script computes LD scores for the custom annotations using ldsc.py.
# Runs as a SLURM array job with one task per chromosome (1-22).
# Uses GPU partition for faster computation.
# =============================================================================

# Source configuration
SCRIPT_DIR="/ocean/projects/bio250020p/kbenjamin/projects/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

# Current chromosome from array task ID
CHR=${SLURM_ARRAY_TASK_ID}
log_message "Processing chromosome: $CHR"

# Validate resources
if ! validate_resources; then
    log_message "ERROR: Resource validation failed. Exiting."
    exit 1
fi

# LDSC script
LDSC_SCRIPT="${LDSC_DIR}/ldsc.py"

# Resource paths
BIM_DIR="${RESOURCE_DIR}/1000G_Phase3_plinkfiles"
HAPMAP3_SNPS="${RESOURCE_DIR}/w_hm3.snplist"

# Process each brain region and heritability status
for REGION in "${BRAIN_REGIONS[@]}"; do
    log_message "Processing region: $REGION"

    for STATUS in "${HERITABILITY[@]}"; do
        log_message "  Processing status: $STATUS"

        # Annotation and output directory
        OUT_DIR="./custom_ldscores/${REGION}/${STATUS}"

        # Input annotation file
        ANNOT_FILE="${OUT_DIR}/${REGION}_${STATUS}.${CHR}.annot.gz"

        # Output file prefix
        OUT_PREFIX="${OUT_DIR}/${REGION}_${STATUS}.${CHR}"

        # Check if annotation file exists
        if [[ ! -f "$ANNOT_FILE" ]]; then
            log_message "    WARNING: Annotation file not found: $ANNOT_FILE"
            continue
        fi

        # Skip if LD scores already exist
        if [[ -f "${OUT_PREFIX}.l2.ldscore.gz" ]]; then
            log_message "    LD scores already exist, skipping..."
            continue
        fi

        log_message "    Computing LD scores..."
        python "$LDSC_SCRIPT" \
            --l2 \
            --bfile "${BIM_DIR}/1000G.EUR.QC.${CHR}" \
            --ld-wind-cm 1 \
            --annot "$ANNOT_FILE" \
            --thin-annot \
            --out "$OUT_PREFIX" \
            --print-snps "$HAPMAP3_SNPS"
    done
done

log_message "**** Job ends ****"
