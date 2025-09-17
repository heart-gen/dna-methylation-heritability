#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=format_expression
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --output=logs/formatting.%j.log
#SBATCH --time=00:30:00

# Function to echo with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility
log_message "**** Loading modules ****"

module purge
module load htslib/1.16
module list

# Set path variables
log_message "**** Loading mamba environment ****"
ENV_PATH="/projects/p32505/opt/env"

mamba run -p $ENV_PATH/AI_env \
      python ../_h/03.prepare_expression.py \
      ./normalized_methylation.tsv \
      ./vcf_chr_list.txt vmrs \
      -o ./ --bed_file ./feature.bed \
      --skip_sample_lookup \
      --sample_id_list ./sample_brnum.txt
      
if [ $? -ne 0 ]; then
    log_message "Error: Python script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
