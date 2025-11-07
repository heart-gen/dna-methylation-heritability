#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=02:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=16G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-24
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=cpg_stats  # Job name
#SBATCH --output=/dev/null      # Standard output log
#SBATCH --error=/dev/null       # Standard error log

# Map chr names 
if [ "$SLURM_ARRAY_TASK_ID" -le 22 ]; then
    CHR="$SLURM_ARRAY_TASK_ID"
elif [ "$SLURM_ARRAY_TASK_ID" -eq 23 ]; then
    CHR="X"
elif [ "$SLURM_ARRAY_TASK_ID" -eq 24 ]; then
    CHR="Y"
else
    echo "Invalid SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Create log directories for each chr
LOG_DIR="logs/chr_${CHR}"
mkdir -p "$LOG_DIR"

# Redirect output and error logs to chr-specific log files
exec > >(tee -a "$LOG_DIR/cpg_stats_${CHR}_out.log")
exec 2> >(tee -a "$LOG_DIR/cpg_stats_${CHR}_err.log" >&2)

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
ENV_PATH="/projects/p32505/opt/envs"

echo "Working on: Chromosome "$CHR

## Activate conda environment
#source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda run -p $ENV_PATH/epigenomics Rscript ../_h/01.get_cpg_stats.R $CHR
if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
