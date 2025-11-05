#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --job-name=munge_sumstats  # Job name
#SBATCH --output=logs/01.output_munge_sumstats.log  # Standard output log
#SBATCH --error=logs/01.error_munge_sumstats.log    # Standard error log

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

SCRIPT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/ldsc
GWAS_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/gwas
HAPMAP3_SNPS=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/hm3_no_MHC.list.txt
OUT_DIR=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/munged_sumstats
mkdir -p $OUT_DIR

python $SCRIPT_DIR/munge_sumstats.py \
    --sumstats $GWAS_DIR/pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv.gz \
    --snp ID \
    --a1 EA \
    --a2 NEA \
    --p PVAL \
    --N-col NCON \
    --signed-sumstats BETA,0 \
    --out $OUT_DIR/MDD_GWAS_munged

zcat "$OUT_DIR/MDD_GWAS_munged.sumstats.gz" | \
awk 'NR==FNR {snps[$1]; next} FNR==1 || $1 in snps' "$HAPMAP3_SNPS" - | \
gzip > "$OUT_DIR/MDD_GWAS_munged_filtered.sumstats.gz"

log_message "**** Job ends ****"