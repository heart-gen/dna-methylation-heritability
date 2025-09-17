#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=prep_covs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --output=summary.%j.log
#SBATCH --time=00:10:00

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
ENV_PATH="/projects/p32505/opt/env"

mamba run -p ${ENV_PATH}/r_env Rscript ../_h/04.generate_covs.R

if [ $? -ne 0 ]; then
    log_message "Error: Rscript execution failed"
    exit 1
fi

log_message "**** Job ends ****"
