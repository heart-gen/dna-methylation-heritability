#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=1G                # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=vmr_test  # Job name
#SBATCH --output=vmr_%j_out.log  # Standard output log
#SBATCH --error=vmr_%j_err.log    # Standard error log

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

## List current modules for reproducibility

module purge
module load R/4.3.0
module list 

# Set path variables
VENV_PATH="/projects/p32505/opt/env/R_env/"

# Activate virtual environment
log_message "**** Activate virtual environment ****"
source "${VENV_PATH}/bin/activate" 

if [ $? -ne 0 ]; then
    log_message "Error: Failed to activate virtual environment"
    exit 1
fi

# Run main script
Rscript --verbose ../_h/01.sd_test.R

if [ $? -ne 0 ]; then
    log_message "Error: R script execution failed"
    exit 1
fi

deactivate