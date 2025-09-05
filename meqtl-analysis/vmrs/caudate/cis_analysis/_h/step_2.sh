#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --job-name=permutation_eqtl
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sierramannion2028@u.northwestern.edu ## Update!
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --output=logs/permutation.%j.log
#SBATCH --time=2:00:00

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
module load cuda/11.6.2-gcc-12.3.0
module list

# Set path variables
log_message "**** Loading mamba environment ****" ## Update from mamba to conda!
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
eval "$(mamba shell hook --shell bash)"

ENV_PATH="/projects/p32505/opt/env/eQTL_env"

mamba run -p ${ENV_PATH} python ../_h/02.permutation_analysis.py

if [ $? -ne 0 ]; then
    log_message "Error: mamba or script execution failed"; exit 1
fi

log_message "**** Job ends ****"
