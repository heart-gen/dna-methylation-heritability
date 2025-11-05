#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --job-name=munge_sumstats  # Job name
#SBATCH --output=logs/output_munge_sumstats.log  # Standard output log
#SBATCH --error=logs/error_munge_sumstats.log    # Standard error log

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

mkdir -p munged_sumstats

python /projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/ldsc/munge_sumstats.py \
    --sumstats /projects/b1213/resources/public_data/gwas/alz/AD_sumstats_Jansenetal_2019sept.txt.gz \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --N-col Nsum \
    --signed-sumstats BETA,0 \
    --out munged_sumstats/alz

log_message "**** Job ends ****"