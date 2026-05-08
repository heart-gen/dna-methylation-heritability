#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=normal
#SBATCH --time=08:00:00
#SBATCH --mem=10gb
#SBATCH --array=1-22
#SBATCH --job-name=compute_ldscores
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/compute_ldscores.%A_%a.log

# =============================================================================
# Step 3: Compute LD Scores
# =============================================================================
# This script computes LD scores for the custom annotations using ldsc.py.
# Runs as a SLURM array job with one task per chromosome (1-22).
# =============================================================================

# Source configuration
SCRIPT_DIR="/projects/b1213/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/all_individuals/tissue_comparison/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

module purge
module list

log_message "**** Loading conda environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

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

    # Annotation and output directory
    OUT_DIR="./custom_ldscores/${REGION}/EA"

    # Input annotation file
    ANNOT_FILE="${OUT_DIR}/${REGION}.${CHR}.annot.gz"

    # Output file prefix
    OUT_PREFIX="${OUT_DIR}/${REGION}.${CHR}"

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
    if ! python "$LDSC_WRAPPER" "$LDSC_DIR" ldsc.py \
        --l2 \
        --bfile "${BIM_DIR}/1000G.EUR.QC.${CHR}" \
        --ld-wind-cm 1 \
        --annot "$ANNOT_FILE" \
        --thin-annot \
        --out "$OUT_PREFIX" \
        --print-snps "$HAPMAP3_SNPS"; then
        log_message "    ERROR: LD score computation failed for $REGION chr$CHR"
        continue
    fi
done

conda deactivate
log_message "**** Job ends ****"
