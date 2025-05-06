#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=10G                # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=plot_enrichment # Job name
#SBATCH --output=logs/plot_enrich.%j.log # Standard output log

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
ENV_PATH="/projects/p32505/opt/env"

log_message "Plotting GO and KEGG enrichment for each brain region"

## Activate conda environment
conda run -p $ENV_PATH/R_env Rscript ../_h/02.plot_enrichment.R

if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

