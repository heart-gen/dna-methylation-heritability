#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=5G                # Memory limit
#SBATCH --job-name=liftOver  # Job name
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

INPUT_FILE=$1

zcat $INPUT_FILE | \
awk 'NR>1 {print $1"\t"$2-1"\t"$3"\t"$4}' \
> gwas.bed

# Runs liftOver
conda activate /projects/p32505/opt/envs/genomics
liftOver gwas.bed ../hg19ToHg38.over.chain.gz lifted.bed gwas.bed
conda deactivate

awk '{print $1"\t"$2}' lifted.bed > new_positions.txt

zcat $INPUT_FILE | \
awk 'BEGIN{OFS="\t"}
NR==FNR {pos[$1]=$2; next}
NR==1 {print; next}
{
    if($1 in pos) $4=pos[$1];
    print
}' new_positions.txt - | \
gzip > gwas_lifted.txt.gz

log_message "**** Job ends ****"