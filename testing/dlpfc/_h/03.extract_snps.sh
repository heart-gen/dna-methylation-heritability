#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics-gpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=16G
#SBATCH --array=0-57
#SBATCH --job-name=extract_snps_chunks
#SBATCH --output=logs/chr_21/extract_snps_chunks_out.log
#SBATCH --error=logs/chr_21/extract_snps_chunks_err.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}
log_message "**** Job starts ****"

# --- Static params ---
CHR="21"
REGION_DIR="./chunked_cpg/chr_${CHR}"
CHROM_SIZES="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"
BASE_NAME="/projects/b1213/users/alexis/projects/dna-methylation-heritability/inputs/genotypes/TOPMed_LIBD.AA"
SAMPLE_LIST="./samples.txt"
OUTPUT_DIR="plink_outputs/chr_${CHR}"
mkdir -p "$OUTPUT_DIR"

# --- Load GNU parallel ---
module load parallel/20160922

# --- Derive chunk file ---
CHUNK_FILE=$(printf "$REGION_DIR/chunk_%03d.bed" "$SLURM_ARRAY_TASK_ID")
if [[ ! -f "$CHUNK_FILE" ]]; then
    echo "Missing region file: $CHUNK_FILE"
    exit 1
fi

# --- Parallel PLINK extraction ---
parallel -j ${SLURM_CPUS_PER_TASK} --colsep '\t' "
    plink2 \
      --pfile $BASE_NAME \
      --chr $CHR \
      --extract range <(echo {1} '\t' {2} '\t' {3}) \
      --make-bed \
      --keep $SAMPLE_LIST \
      --no-parents \
      --no-sex \
      --no-pheno \
      --out $OUTPUT_DIR/chr${CHR}_region_{#}_snps"

log_message "**** Job ends ****"
