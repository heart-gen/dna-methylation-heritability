#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics-gpu
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --array=1-22
#SBATCH --job-name=compute_ldscores
#SBATCH --gres=gpu:a100:1
#SBATCH --output=logs/03.output_compute_ldscores.%a.log
#SBATCH --error=logs/03.error_compute_ldscores.%a.log

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

SCRIPT=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/ldsc/ldsc.py
BIM_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/GRCh38/plink_files
HAPMAP3_SNPS=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/hm3_no_MHC.list.txt

for REGION in "${BRAIN_REGIONS[@]}"; do
	for STATUS in "${HERITABILITY[@]}"; do

    	OUT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/custom_ldscores/${REGION}/${STATUS}

    	python $SCRIPT \
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