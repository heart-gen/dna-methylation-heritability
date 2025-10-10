#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:1
#SBATCH --job-name=cis_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.bennett@northwestern.edu ## Update this
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --output=logs/nominal-analysis.%j.log
#SBATCH --time=00:20:00

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
log_message "**** Loading conda environment ****"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
python ../_h/01.eqtl_tensorqtl.py

#ENV_PATH="/projects/p32505/opt/env/eQTL_env"
#conda run -p ${ENV_PATH} python ../_h/01.eqtl_tensorqtl.py

if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

conda deactivate

log_message "**** Job ends ****"
