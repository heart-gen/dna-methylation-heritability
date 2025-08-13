#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics-gpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=40G
#SBATCH --job-name=extract_snps
#SBATCH --output=logs/extract_snps_out.log
#SBATCH --error=logs/extract_snps_err.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

# --- Static params ---
CHROM_SIZES="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"
BASE_NAME="/projects/b1213/users/alexis/projects/dna-methylation-heritability/inputs/genotypes/TOPMed_LIBD.AA"
SAMPLE_LIST="./samples.txt"

# --- Load GNU parallel ---
module load parallel/20160922

# --- Map global SLURM_ARRAY_TASK_ID to chromosome & chunk ---
# Mapping file must exist
MAP_FILE="chunk_map.tsv"
if [[ ! -f "$MAP_FILE" ]]; then
    echo "Mapping file $MAP_FILE not found!" >&2
    exit 1
fi

# Get line corresponding to current task
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MAP_FILE")
CHR=$(echo "$line" | cut -f1)
CHUNK_FILE=$(echo "$line" | cut -f2)

if [[ -z "$CHR" || -z "$CHUNK_FILE" ]]; then
    echo "Invalid mapping for array ID $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

# --- Output & log dirs ---
OUTPUT_DIR="plink_outputs/chr_${CHR}"
mkdir -p "$OUTPUT_DIR"

# --- Parallel PLINK extraction ---
parallel -j ${SLURM_CPUS_PER_TASK} --colsep '\t' "
    tmpfile=\$(mktemp)
    echo -e \"{1}\t{2}\t{3}\" > \"\$tmpfile\"
    plink2 \
      --pfile $BASE_NAME \
      --chr $CHR \
      --extract range \$tmpfile \
      --make-bed \
      --keep $SAMPLE_LIST \
      --no-parents \
      --no-sex \
      --no-pheno \
      --out $OUTPUT_DIR/chr${CHR}_region_{#}_snps
    rm -f \"\$tmpfile\"" :::: "$CHUNK_FILE"

log_message "**** Job ends ****"
