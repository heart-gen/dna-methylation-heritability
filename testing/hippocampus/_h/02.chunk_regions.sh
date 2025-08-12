#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --job-name=chunk_bed
#SBATCH --output=logs/chunk_regions_out.log
#SBATCH --error=logs/chunk_regions_err.log

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

# Parameters
CHUNK_SIZE=5000
BED_FILENAME="formatted_cpg.bed"
INPUT_BASE_DIR="./cpg"
OUTPUT_BASE_DIR="./chunked_cpg"

mkdir -p "$OUTPUT_BASE_DIR"

# Loop through chromosomes 1–22
for CHR in {1..22}; do
    CHR_DIR="${INPUT_BASE_DIR}/chr_${CHR}"
    INPUT_BED="${CHR_DIR}/${BED_FILENAME}"
    CHR_OUTPUT_DIR="${OUTPUT_BASE_DIR}/chr_${CHR}"

    mkdir -p "$CHR_OUTPUT_DIR"

    # Split into chunks with padded suffix
    split -l "$CHUNK_SIZE" -d -a 3 --additional-suffix=".bed" \
        "$INPUT_BED" "${CHR_OUTPUT_DIR}/chunk_"

    # Remove the filtered intermediate file
    rm "$INPUT_BED"
done

log_message "**** Job ends ****"