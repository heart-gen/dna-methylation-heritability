#!/bin/bash
#SBATCH --account=b1213        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=1G                # Memory limit
#SBATCH --job-name=gwas  # Job name
#SBATCH --output=logs/output_%j.log  # Standard output log
#SBATCH --error=logs/error_%j.log    # Standard error log

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

GWAS_FILE=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/liftOver/_m/AFR_lifted.sumstats.txt
REGIONS=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/fine-mapping/_m/cal_ld/chr_21

for dir in $REGIONS/*; do
    region=$(basename "$dir")
    start=${region%_*}
    end=${region#*_}

    awk -v s=$start -v e=$end '
    NR==1 {print; next}
    {
        split($2,a,":")
        if (a[2]>=s && a[2]<=e) print
    }' $GWAS_FILE > "$dir/gwas_$start_$end.txt"
done

log_message "**** Job ends ****"