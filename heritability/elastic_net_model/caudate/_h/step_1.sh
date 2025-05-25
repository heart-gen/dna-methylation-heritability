#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=elastic_h2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/elastic_h2_%A_%a.out
#SBATCH --error=logs/elastic_h2_%A_%a.err
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** BRIDGES info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## List current modules for reproducibility

module purge
module load anaconda3/2024.10-1
module list 

# Activate project conda env
ENV_PATH="/ocean/projects/bio250020p/shared/opt/env"

# Define region if applicable
export region="caudate"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

log_message "**** Run elastic net ****"
conda run -p $ENV_PATH/R_env Rscript ../_h/01.elastic-net.R

log_message "**** Job ends ****"
