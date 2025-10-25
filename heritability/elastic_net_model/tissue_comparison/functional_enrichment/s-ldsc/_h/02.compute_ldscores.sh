#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --array=1-22            # Array job for chromosomes 1-22
#SBATCH --job-name=compute_ldscores  # Job name
#SBATCH --output=logs/compute_ldscores/output_%j.log  # Standard output log
#SBATCH --error=logs/compute_ldscores/error_%j.log    # Standard error log

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

output_dir="custom_ldscores"
mkdir -p ${output_dir}

python ldsc.py \
	--l2 \
	--bfile /projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/plink_files/chr_${SLURM_ARRAY_TASK_ID} \
	--ld-wind-cm 1 \
	--annot /projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/annotations/chr_${SLURM_ARRAY_TASK_ID}.annot.gz \
	--thin-annot \
	--out $output_dir \
	--print-snps /projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/w_hm3.snplist.gz

log_message "**** Job ends ****"