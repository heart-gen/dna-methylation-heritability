#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=00:10:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=5G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=summarize  # Job name
#SBATCH --output=summarize.%j.log      # Standard output log

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

## Edit with your job command
echo "**** Summarize sample data ****"

## Activate conda environment
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/epigenomics
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
python ../_h/02.summarize.py
conda deactivate

if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
