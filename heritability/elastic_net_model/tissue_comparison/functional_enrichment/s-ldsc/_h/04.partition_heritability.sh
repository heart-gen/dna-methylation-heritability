#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name=partition_heritability
#SBATCH --output=logs/04.output_partition_heritability.log
#SBATCH --error=logs/04.error_partition_heritability.log

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

BRAIN_REGIONS=("caudate" "dlpfc" "hippocampus")
HERITABILITY=("heritable" "non_heritable" "low_prediction")
DISEASES=("ad" "mdd" "scz")

SCRIPT=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/ldsc/ldsc.py
SUMSTATS_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/sumstats
BASELINE_LD_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/GRCh38/baselineLD_v2.2
WEIGHTS_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/GRCh38/weights
FRQ_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/GRCh38/plink_files

for DISEASE in "${DISEASES[@]}"; do
	echo "Processing disease: ${DISEASE}"

	for REGION in "${BRAIN_REGIONS[@]}"; do

		echo "Processing region: ${REGION}"

		for STATUS in "${HERITABILITY[@]}"; do

			echo "Processing status: ${STATUS}"

			CUSTOM_LD_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/custom_ldscores/${REGION}/${STATUS}
			OUT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/results/${DISEASE}/${REGION}/${STATUS}

			mkdir -p $OUT_DIR

			for CHR in {1..22}; do

			echo "Processing chromosome: ${CHR}"

				python $SCRIPT \
		    		--h2 $SUMSTATS_DIR/${DISEASE}.sumstats.gz \
		    		--ref-ld $BASELINE_LD_DIR/baselineLD.${CHR} \
		    		--ref-ld $CUSTOM_LD_DIR/${REGION}_${STATUS}.${CHR} \
		    		--w-ld $WEIGHTS_DIR/weights.hm3_noMHC.${CHR} \
					--frqfile $FRQ_DIR/1000G.EUR.hg38.${CHR} \
					--overlap-annot \
					--thin-annot \
		    		--out $OUT_DIR/${DISEASE}_${REGION}_${STATUS}.${CHR}
			done
		done
	done
done

log_message "**** Job ends ****"
