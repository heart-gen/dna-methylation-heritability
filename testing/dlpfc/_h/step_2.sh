#!/bin/bash
#SBATCH --job-name=reformat_cpgs
#SBATCH --output=logs/02.reformat_cpgs_%j.out
#SBATCH --error=logs/02.reformat_cpgs_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --account=p32505
#SBATCH --mail-user=elisajohnson2027@u.northwestern

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

## List current modules for reproducibility

module purge
module list

# Set path variables
ENV_PATH="/projects/p32505/opt/env"

## Activate conda environment
#source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda run -p $ENV_PATH/r_env Rscript ../_h/02.reformat_cpgs.R
if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

log_message "**** Job ends ****"