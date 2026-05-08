#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=music_cellprop
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH --output=logs/music_cellprop.%A_%a.log
#SBATCH --array=0-2
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
echo "Node name: ${SLURM_NODENAME:-N/A}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility
log_message "**** Loading modules ****"

module purge
module list

# Set path variables
log_message "**** Loading conda R environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

REGIONS=(caudate dlpfc hippocampus)
REGION="${REGIONS[SLURM_ARRAY_TASK_ID]}" || { log_message "Error: Unknown array index ${SLURM_ARRAY_TASK_ID}"; exit 1; }

# Run main script
log_message "**** Run Rscript ****"
Rscript ../_h/01.estimate_cell_prop.R "${REGION}"

if [ $? -ne 0 ]; then
    log_message "Error: Rscript execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
