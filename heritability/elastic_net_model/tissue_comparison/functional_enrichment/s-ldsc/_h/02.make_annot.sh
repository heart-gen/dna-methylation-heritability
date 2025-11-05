#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --job-name=make_annot  # Job name
#SBATCH --output=logs/02.output_make_annot.log  # Standard output log
#SBATCH --error=logs/02.error_make_annot_.log    # Standard error log

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

module load bedtools/2.30.0

BRAIN_REGIONS=("caudate" "dlpfc" "hippocampus")
HERITABILITY=("heritable" "non_heritable" "low_prediction")

SCRIPT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/ldsc
BIM_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/GRCh38/plink_files

for REGION in "${BRAIN_REGIONS[@]}"; do

  echo "Processing region: $REGION"

  for STATUS in "${HERITABILITY[@]}"; do
    echo "Processing status: $STATUS"

    BED_FILE=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/vmr/${REGION}/${STATUS}.bed
    OUT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/custom_ldscores/${REGION}/${STATUS}
    mkdir -p $OUT_DIR

    for CHR in {1..22}; do
      echo "Processing chromosome $CHR..."

      BIM_FILE=${BIM_DIR}/1000G.EUR.hg38.${CHR}.bim

      python $SCRIPT_DIR/make_annot.py \
        --bed-file ${BED_FILE} \
        --bimfile ${BIM_FILE} \
        --annot-file $OUT_DIR/${REGION}_${STATUS}.${CHR}.annot.gz \
        --windowsize 500000 \

    done
  done
done

log_message "**** Job ends ****"