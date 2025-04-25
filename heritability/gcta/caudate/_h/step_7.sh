#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=5G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=combine_results      # Job name
#SBATCH --output=/logs/combine_results_$j_out.log      # Standard output log
#SBATCH --error=/logs/combine_results_$j_err.log       # Standard error log

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

## Edit with your job command
echo "Aggregating GREML output for all analyzed regions"

## Activate conda environment
conda run -p $ENV_PATH/AI_env python ../_h/07.combine_results.py

if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
