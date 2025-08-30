#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --array=0-8118
#SBATCH --job-name=extract_snps
#SBATCH --output=logs/extract_snps_%A_%a.out
#SBATCH --error=logs/extract_snps_%A_%a.err

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}
log_message "**** Job starts ****"

## List current modules for reproducibility

module purge
module load plink/2.0-alpha-3.3
module list 

# --- Derive chunk file ---
CHUNK_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" chunk_file_list.txt)

if [[ ! -f "$CHUNK_FILE" ]]; then
    echo "Missing region file: $CHUNK_FILE"
    exit 1
fi

# Extract CHR from path: cpg/chr_1/chunked_cpg/chunk_001.bed
CHR=$(echo "$CHUNK_FILE" | sed -E 's|.*chr_([0-9XY]+).*|\1|')

# --- Static params ---
BASE_NAME="../../../inputs/genotypes/TOPMed_LIBD.AA"
SAMPLE_LIST="../../../heritability/caudate/_m/samples.txt"
OUTPUT_DIR="./plink_outputs/chr_${CHR}"
mkdir -p "$OUTPUT_DIR"

log_message "Processing chunk: $CHUNK_FILE"

while IFS=$'\t' read -r chr start end label; do
    OUT_FILE="${OUTPUT_DIR}/cpg_TOPMed_LIBD.AA.${start}_${end}"

    log_message "Extracting SNPs for $chr:$start-$end → $OUT_FILE"

    plink2 \
        --pfile "$BASE_NAME" \
        --chr "$CHR" \
        --from-bp "$start" \
        --to-bp "$end" \
        --make-bed \
        --keep "$SAMPLE_LIST" \
        --no-parents \
        --no-sex \
        --no-pheno \
        --out "$OUT_FILE"
done < "$CHUNK_FILE"

log_message "**** Job ends ****"