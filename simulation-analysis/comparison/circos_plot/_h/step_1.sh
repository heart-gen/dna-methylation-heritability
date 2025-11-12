#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=plot_circos
#SBATCH --mail-type=FAIL ## If you want to have it email you for any reason
#SBATCH --mail-user=alexis.bennett@northwestern.edu ## replace with your email
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --output=logs/plot_circos.%j.log
#SBATCH --time=00:05:00

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
module list

# Set path variables
log_message "**** Loading mamba environment ****"
ENV_PATH="/projects/p32505/opt/envs"

mamba run -p $ENV_PATH/epigenomics Rscript ../_h/01.plot_circos.R ## change to conda if mamba errors
if [ $? -ne 0 ]; then
    log_message "Error: Mamba or script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
