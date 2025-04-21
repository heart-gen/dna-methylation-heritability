#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=15G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-12078%300
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=cal_vmr  # Job name
#SBATCH --output=logs/cal_vmr_%j_out.log  # Standard output log
#SBATCH --error=logs/cal_vmr_%j_err.log    # Standard error log

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
module load python/3.10.1
module list 

## Edit with your job command
REGION_LIST="./vmr_list.txt"

# Get the current sample name from the sample list
REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $REGION_LIST)
CHR=$(echo "$REGION" | awk '{print $1}')
START=$(echo "$REGION" | awk '{print $2}')
END=$(echo "$REGION" | awk '{print $3}')

# Set path variables
ENV_PATH="/projects/p32505/opt/env"

echo "Working on: Chromosome "$CHR:$START-$END 

# Run main script
python ../_h/06.summary.py

if [ $? -ne 0 ]; then
    log_message "Error: Python script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
