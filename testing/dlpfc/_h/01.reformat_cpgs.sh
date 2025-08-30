#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --array=1-22
#SBATCH --job-name=reformat_cpgs
#SBATCH --output=logs/reformat_cpgs_%A_%a.out
#SBATCH --error=logs/reformat_cpgs_%A_%a.err

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

CHR="${SLURM_ARRAY_TASK_ID}"
CPG_BASE_DIR="../../../heritability/dlpfc/_m/cpg/chr_${CHR}"
INPUT_FILE="${CPG_BASE_DIR}/cpg_pos.txt"
OUTPUT_DIR="./cpg/chr_$CHR"
mkdir -p "$OUTPUT_DIR"
OUTPUT_FILE="${OUTPUT_DIR}/formatted_cpg.bed"

#input_file_name="cpg_pos.txt"
chrom_sizes_file="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"

if [[ ! -f "$INPUT_FILE" ]]; then
    log_message "[SKIP] chr_$CHR — $INPUT_FILE not found."
    exit 0
fi

if [[ ! -f "$chrom_sizes_file" ]]; then
    log_message "Missing chromosome size file: $chrom_sizes_file"
    exit 1
fi

# Load chromosome sizes, stripping "chr" prefix for array keys
declare -A chrom_sizes
while read -r chrom size; do
    chrom_no="${chrom#chr}"   # remove 'chr' prefix
    chrom_sizes["$chrom_no"]=$size
done < "$chrom_sizes_file"

chrom_size=${chrom_sizes["$CHR"]}
if [[ -z "$chrom_size" ]]; then
    log_message "[SKIP] Chromosome size not found for $CHR."
    exit 0
fi

log_message "[INFO] Processing chromosome $CHR ($INPUT_FILE)..."

row=1
site_count=0
skipped_count=0

> "$OUTPUT_FILE"

while IFS= read -r X; do
    bed_start=$((X - 20000))      # no -1, keep 1-based
    bed_end=$((X + 20000))        # inclusive end

    # Skip if out of bounds
    if (( bed_start <= 0 || bed_end >= chrom_size )); then
        skipped_count=$((skipped_count + 1))
        continue
    fi

    echo -e "chr${CHR}\t${bed_start}\t${bed_end}\tCpG_${row}" >> "$OUTPUT_FILE"
    row=$((row + 1))
    site_count=$((site_count + 1))
done < <(tail -n +3 "$INPUT_FILE")

log_message "[DONE] Wrote $site_count CpG sites to $OUTPUT_FILE (skipped $skipped_count out-of-bounds)"

log_message "**** Job ends ****"
