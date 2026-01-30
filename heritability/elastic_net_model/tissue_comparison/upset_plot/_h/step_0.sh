#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=5
#SBATCH --job-name=region_heritability
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --output=logs/region_heritability.%j.log

# =============================================================================
# Step 0: Separate VMRs by heritability status
# =============================================================================
# This script processes elastic-net summary files for each brain region and
# separates VMRs into heritable and non-heritable categories
# =============================================================================

# Source configuration
source "../_h/config.sh"

log_message "**** Job starts ****"
print_job_info

# Create logs directory if it doesn't exist
mkdir -p logs

# Python script for region heritability separation
PYTHON_SCRIPT="../_h/00.region_heritability.py"

# Chain file for liftover (hg38 to hg19)
CHAIN_FILE="../../../../../../../inputs/supportfiles/_m/hg38ToHg19.over.chain"

# Process each brain region
for REGION in "${BRAIN_REGIONS[@]}"; do
    log_message "Processing region: $REGION"

    # Input file path
    INPUT_FILE="../../../../../${REGION}/_m/${REGION}_summary_elastic-net.tsv"

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

log_message "**** Job ends ****"
