#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=00:10:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=8G                # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=upset_plot # Job name
#SBATCH --output=logs/upset.%j.log # Standard output log

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
module list

# Set path variables
ENV_PATH="/projects/p32505/opt/env/intervene_env"
OVERLAP=(0.25 0.5 0.75)

for f in "${OVERLAP[@]}"; do
    OUTDIR="./f_${f}"
    mkdir -p "$OUTDIR"

    log_message "Analyzing brain region overlap using Intervene for f=$f"
    
    ## Activate environment
    conda run -p "$ENV_PATH" intervene upset \
      -i ./*.bed \
      --bedtools-options f=$f \
      --save-overlaps \
      --output "$OUTDIR"

    if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
    fi
done

for F in "${OVERLAP[@]}"; do
    OUTDIR="./F_${F}"
    mkdir -p "$OUTDIR"

    log_message "Analyzing brain region overlap using Intervene for F=$F"
    
    ## Activate environment
    conda run -p "$ENV_PATH" intervene upset \
      -i ./*.bed \
      --bedtools-options F=$F \
      --save-overlaps \
      --output "$OUTDIR"

    if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
    fi
done

log_message "**** Job ends ****"