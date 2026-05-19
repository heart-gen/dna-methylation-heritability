#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=region_heritability
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/region_heritability.%j.log

# =============================================================================
# Step 0: Separate VMRs by heritability status and liftover to hg19
# =============================================================================
# This script processes elastic-net summary files for each brain region,
# separates VMRs into heritable and non-heritable categories, and performs
# liftover from hg38 to hg19 coordinates.
# =============================================================================

# Source configuration
SCRIPT_DIR="/projects/b1213/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/BA_only/tissue_comparison/sensitivity_analysis/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

module purge
module list

log_message "**** Loading conda environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

# Python script for region heritability separation
PYTHON_SCRIPT="${SCRIPT_DIR}/region_heritability.py"

# Process each brain region
for REGION in "${BRAIN_REGIONS[@]}"; do
    log_message "Processing region: $REGION"

    # Input file path
    INPUT_FILE="$(get_input_file "$REGION")"

    # Output directory
    OUTPUT_DIR="./vmr/${REGION}"
    mkdir -p "$OUTPUT_DIR"

    # Verify input file exists
    if [[ ! -f "$INPUT_FILE" ]]; then
        log_message "ERROR: Input file not found: $INPUT_FILE"
        continue
    fi

    log_message "  Input: $INPUT_FILE"
    log_message "  Output: $OUTPUT_DIR"
    log_message "  Chain file: $CHAIN_FILE"

    # Run separation script
    python "$PYTHON_SCRIPT" \
        --input_file "$INPUT_FILE" \
        --output_dir "$OUTPUT_DIR" \
        --chain_file "$CHAIN_FILE"
done

conda deactivate
log_message "**** Job ends ****"
