#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --job-name=chunk_bed
#SBATCH --array=1-22
#SBATCH --output=logs/chunk_regions_%A_%a.out
#SBATCH --error=logs/chunk_regions_%A_%a.err

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
CHR="${SLURM_ARRAY_TASK_ID}"
CHUNK_SIZE=3000
BED_FILENAME="formatted_cpg.bed"
BASE_DIR="./cpg/chr_$CHR"
INPUT_BED="${BASE_DIR}/${BED_FILENAME}"
OUTPUT_BASE_DIR="${BASE_DIR}/chunked_cpg"

mkdir -p "$OUTPUT_BASE_DIR"

log_message "Splitting chr$CHR bed into chunks of $CHUNK_SIZE lines:"

# Split into chunks with padded suffix
split -l "$CHUNK_SIZE" -d -a 3 --additional-suffix=".bed" \
    "$INPUT_BED" "${OUTPUT_BASE_DIR}/chunk_"

log_message "**** Job ends ****"