#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --array=1-22            # Array job for chromosomes 1-22
#SBATCH --job-name=make_annot  # Job name
#SBATCH --output=logs/make_annot/output_%j.log  # Standard output log
#SBATCH --error=logs/make_annot/error_%j.log    # Standard error log

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

module load bedtools/2.30.0

output_dir="custom_annotations"
mkdir -p ${output_dir}

echo Processing chromosome ${SLURM_ARRAY_TASK_ID}

python /projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/ldsc/make_annot.py \
    --bed-file /projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/plink_files/chr_${SLURM_ARRAY_TASK_ID}.bed \
    --bim /projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/_m/plink_files/chr_${SLURM_ARRAY_TASK_ID}.bim \
    --annot-file annotations/chr_${SLURM_ARRAY_TASK_ID}.annot.gz \
    --windowsize 500000

echo Annotation for chromosome ${SLURM_ARRAY_TASK_ID} created.

log_message "**** Job ends ****"