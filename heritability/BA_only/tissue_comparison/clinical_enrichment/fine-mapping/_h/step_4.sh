#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=20gb
#SBATCH --job-name=compute_ld_matrix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/compute_ld_matrix.%j.log
#SBATCH --array=1-22

# Function to echo with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

module load plink/1.9

PLINK_DIR=/projects/b1213/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/fine-mapping/_m/plink/hippocampus/chr_$SLURM_ARRAY_TASK_ID
OUT_DIR=./ld_matrices/hippocampus/chr_$SLURM_ARRAY_TASK_ID
CHR=$SLURM_ARRAY_TASK_ID

mkdir -p ${OUT_DIR}

echo "Input directory: ${PLINK_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "Chromosome: ${CHR}"
echo ""

for BED in ${PLINK_DIR}/*.bed; do

    PREFIX=$(basename ${BED} .bed)
    FULL_PREFIX=${PLINK_DIR}/${PREFIX}
    OUT_PREFIX=${OUT_DIR}/${PREFIX}

    # Extract START and END from filename
    # Assumes pattern: something.START_END_something
    COORD_PART=$(echo ${PREFIX} | grep -o '[0-9]\+_[0-9]\+')
    START=$(echo ${COORD_PART} | cut -d'_' -f1)
    END=$(echo ${COORD_PART} | cut -d'_' -f2)

    if [[ -z "$START" || -z "$END" ]]; then
        echo "Could not extract coordinates from ${PREFIX}, skipping."
        continue
    fi

    # Skip if already computed
    if [[ -f "${OUT_PREFIX}.ld.gz" ]]; then
        echo "Skipping ${PREFIX} (already computed)"
        continue
    fi

    echo "Processing ${PREFIX}"
    echo "  Region: chr${CHR}:${START}-${END}"

    plink \
        --bfile ${FULL_PREFIX} \
        --chr ${CHR} \
        --r square \
        --out ${OUT_PREFIX}

done

echo ""
echo "All regions complete."

log_message "**** Job ends ****"