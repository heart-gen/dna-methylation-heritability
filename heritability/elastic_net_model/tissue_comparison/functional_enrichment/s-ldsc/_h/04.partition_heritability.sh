#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=30G                # Memory limit
#SBATCH --job-name=partition_heritability  # Job name
#SBATCH --output=logs/04.output_partition_heritability.log  # Standard output log
#SBATCH --error=logs/04.error_partition_heritability.log    # Standard error log

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

SCRIPT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/ldsc
SUMSTATS_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/munged_sumstats
BASELINE_LD_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/GRCh38/baselineLD_v2.2
WEIGHTS_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/GRCh38/weights
FRQ_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/GRCh38/plink_files

for REGION in "${BRAIN_REGIONS[@]}"; do

	for STATUS in "${HERITABILITY[@]}"; do

		CUSTOM_LD_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/custom_ldscores/${REGION}/${STATUS}
		OUT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/results/baseline

		mkdir -p $OUT_DIR

		python $SCRIPT_DIR/ldsc.py \
	    	--h2 $SUMSTATS_DIR/SCZ_GWAS_munged_filtered.sumstats.gz \
	    	--ref-ld-chr $BASELINE_LD_DIR/baselineLD. \
	    	--w-ld-chr $WEIGHTS_DIR/weights.hm3_noMHC. \
			--frqfile-chr $FRQ_DIR/1000G.EUR.hg38. \
			--thin-annot \
			--overlap-annot \
	    	--out $OUT_DIR/SCZ_baseline
	done
done

log_message "**** Job ends ****"
