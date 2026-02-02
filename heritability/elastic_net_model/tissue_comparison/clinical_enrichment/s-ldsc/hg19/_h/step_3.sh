#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --time=08:00:00
#SBATCH --ntasks-per-node=5
#SBATCH --array=1-22
#SBATCH --job-name=compute_ldscores
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --output=logs/compute_ldscores.%A_%a.log

# =============================================================================
# Step 3: Compute LD Scores
# =============================================================================
# This script computes LD scores for the custom annotations using ldsc.py.
# Runs as a SLURM array job with one task per chromosome (1-22).
# =============================================================================

# Source configuration
SCRIPT_DIR="/ocean/projects/bio250020p/kbenjamin/projects/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics

# Current chromosome from array task ID
CHR=${SLURM_ARRAY_TASK_ID}
log_message "Processing chromosome: $CHR"

# Validate resources
if ! validate_resources; then
    log_message "ERROR: Resource validation failed. Exiting."
    exit 1
fi

# LDSC wrapper (patches pandas/Python 3 compatibility issues)
LDSC_WRAPPER="${SCRIPT_DIR}/ldsc_wrapper.py"

# Resource paths
BIM_DIR="${RESOURCE_DIR}/1000G_EUR_Phase3_plink"
HAPMAP3_SNPS="${RESOURCE_DIR}/hm3_no_MHC.list.txt"

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
        python "$LDSC_WRAPPER" "$LDSC_DIR" ldsc.py \
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
