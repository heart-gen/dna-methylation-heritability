#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=10G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-22
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=sd_test  # Job name
##SBATCH --output=logs/sd_%j_out.log  # Standard output log
##SBATCH --error=logs/sd_%j_err.log   # Standard error log

# Create log directories for each chr
LOG_DIR="logs/chr_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$LOG_DIR"

# Redirect output and error logs to chr-specific log files
exec > >(tee -a "$LOG_DIR/extract_vmr_${SLURM_ARRAY_TASK_ID}_out.log")
exec 2> >(tee -a "$LOG_DIR/extract_vmr_${SLURM_ARRAY_TASK_ID}_err.log" >&2)

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

echo "Working on: Chromosome "$SLURM_ARRAY_TASK_ID

## Activate conda environment
#source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
$ENV_PATH/R_env/bin/Rscript ../_h/01.sd_test.R $SLURM_ARRAY_TASK_ID
if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
