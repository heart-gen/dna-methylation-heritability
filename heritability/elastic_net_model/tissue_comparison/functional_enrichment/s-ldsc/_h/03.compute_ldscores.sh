#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=04:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --array=1-22               # Array job for chromosomes 1-22
#SBATCH --job-name=compute_ldscores  # Job name
#SBATCH --output=logs/03.output_compute_ldscores.%a.log  # Standard output log
#SBATCH --error=logs/03.error_compute_ldscores.%a.log    # Standard error log

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
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

BRAIN_REGIONS=("caudate" "dlpfc" "hippocampus")
HERITABILITY=("heritable" "non_heritable" "low_prediction")

SCRIPT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/ldsc
BIM_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/GRCh38/plink_files
HAPMAP3_SNPS=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/hm3_no_MHC.list.txt

for REGION in "${BRAIN_REGIONS[@]}"; do

	for STATUS in "${HERITABILITY[@]}"; do

    	OUT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/custom_ldscores/${REGION}/${STATUS}

    	BIM_FILE=${BIM_DIR}/1000G.EUR.hg38.${SLURM_ARRAY_TASK_ID}.bim

    	python $SCRIPT_DIR/ldsc.py \
			--l2 \
			--bfile $BIM_DIR/1000G.EUR.hg38.${SLURM_ARRAY_TASK_ID} \
			--ld-wind-cm 1 \
			--annot $OUT_DIR/${REGION}_${STATUS}.${SLURM_ARRAY_TASK_ID}.annot.gz \
			--thin-annot \
			--out $OUT_DIR/${REGION}_${STATUS}.${SLURM_ARRAY_TASK_ID} \
			--print-snps $HAPMAP3_SNPS

	done
done

log_message "**** Job ends ****"