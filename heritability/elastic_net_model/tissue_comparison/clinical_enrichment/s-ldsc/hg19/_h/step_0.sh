#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name=region_heritability
#SBATCH --output=logs/00.output_%j.log
#SBATCH --error=logs/00.error_%j.log

# =============================================================================
# Step 0: Separate VMRs by heritability status and liftover to hg19
# =============================================================================
# This script processes elastic-net summary files for each brain region,
# separates VMRs into heritable and non-heritable categories, and performs
# liftover from hg38 to hg19 coordinates.
# =============================================================================

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

# Create logs directory if it doesn't exist
mkdir -p logs

# Python script for region heritability separation
PYTHON_SCRIPT="${SCRIPT_DIR}/00.region_heritability.py"

# Chain file for liftover (hg38 to hg19)
CHAIN_FILE="/ocean/projects/bio250020p/kbenjamin/projects/DNAm-biomarkers-SCZ/inputs/supportfiles/_m/hg38ToHg19.over.chain"

# Process each brain region
for REGION in "${BRAIN_REGIONS[@]}"; do
    log_message "Processing region: $REGION"

    # Input file path
    INPUT_FILE="/ocean/projects/bio250020p/kbenjamin/projects/dna-methylation-heritability/heritability/elastic_net_model/${REGION}/_m/${REGION}_summary_elastic-net.tsv"

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

    # Verify output files were created
    if [[ -f "${OUTPUT_DIR}/heritable_hg19.bed" ]] && [[ -f "${OUTPUT_DIR}/non_heritable_hg19.bed" ]]; then
        log_message "  Successfully created BED files:"
        log_message "    - heritable_hg19.bed: $(wc -l < "${OUTPUT_DIR}/heritable_hg19.bed") regions"
        log_message "    - non_heritable_hg19.bed: $(wc -l < "${OUTPUT_DIR}/non_heritable_hg19.bed") regions"
    else
        log_message "  WARNING: Expected output files not found"
    fi

done

log_message "**** Job ends ****"
