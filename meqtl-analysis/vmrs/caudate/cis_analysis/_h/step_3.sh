#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#$BATCH --gres=gpu:h100:1
#SBATCH --job-name=post_processing
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu ## UPATE THIS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --output=logs/post_processing.%j.log
#SBATCH --time=0:30:00

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
log_message "**** Loading mamba environment ****" ## Might need to change mamba to conda
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
eval "$(mamba shell hook --shell bash)"

ENV_PATH="/projects/p32505/opt/env/eQTL_env"

mamba run -p ${$ENV_PATH} python ../_h/03.post_processing.py

if [ $? -ne 0 ]; then
    log_message "Error: mamba or script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
