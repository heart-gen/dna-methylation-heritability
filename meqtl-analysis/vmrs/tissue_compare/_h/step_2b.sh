#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=prep_model_variables
#SBATCH --mail-type=FAIL ## If you want to have it email you for any reason
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu ## replace with your email
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=150gb
#SBATCH --output=logs/prep_model_variables.%j.log
#SBATCH --time=01:00:00

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
log_message "**** Loading conda environment ****"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/rnaseq

INPUT="general"
Rscript ../_h/02.compute_mash_model_variables.R $INPUT

if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

conda deactivate

log_message "**** Job ends ****"
