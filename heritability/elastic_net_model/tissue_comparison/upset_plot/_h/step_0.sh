#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=00:10:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=5G                # Memory limit
#SBATCH --job-name=partition_heritability
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --output=logs/partition_heritability.%j.log

# =============================================================================
# Step 0: Separate VMRs by heritability status
# =============================================================================
# This script processes elastic-net summary files for each brain region and
# separates VMRs into heritable and non-heritable categories
# =============================================================================

# Log function
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility

module purge
module list

# Set path variables
ENV_PATH="/projects/p32505/opt/envs"

# Output directory
    OUTPUT_DIR="./bed_files"
    mkdir -p "$OUTPUT_DIR"

BRAIN_REGIONS=(caudate hippocampus dlpfc)

# Process each brain region
for REGION in "${BRAIN_REGIONS[@]}"; do
    log_message "Processing region: $REGION"

	BED="${OUTPUT_DIR}/${REGION}.bed"
	if [ ! -f "$BED" ]; then
		log_message "Sorting bed file of all VMRs for ${REGION}."
		cp ../../../../${REGION}/_m/vmr.bed ${OUTPUT_DIR}/${REGION}_unsorted.bed
		sort -k1,1 -k2,2n ${OUTPUT_DIR}/${REGION}_unsorted.bed > ${OUTPUT_DIR}/${REGION}_all.bed
	else
		log_message "Bed file already exists. Skipping."
	fi

	log_message "Partitioning bed files for each heritability class"

    # Input file path
    INPUT_FILE="../../../../elastic_net_model/${REGION}/_m/${REGION}_summary_elastic-net.tsv"

    # Verify input file exists
    if [[ ! -f "$INPUT_FILE" ]]; then
        log_message "ERROR: Input file not found: $INPUT_FILE"
        continue
    fi

    # Run separation script
    conda run -p $ENV_PATH/genomics python ../_h/00.prepare_bed.py \
        --input_file "$INPUT_FILE" \
        --output_dir "$OUTPUT_DIR" \
		--region "$REGION"

	if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
	fi
done

log_message "**** Job ends ****"
