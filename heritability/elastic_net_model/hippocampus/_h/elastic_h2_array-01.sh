#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=elastic_h2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/elastic_h2_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

OFFSET=${OFFSET:-0} # fallback default to 0
task_id=$((OFFSET + SLURM_ARRAY_TASK_ID - 1))
export task_id

echo "**** BRIDGES info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "OFFSET: ${OFFSET}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Computed task_id: ${task_id}"

module purge
module load anaconda3/2024.10-1
module list

ENV_PATH="/ocean/projects/bio250020p/shared/opt/env"
export region="hippocampus"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

log_message "**** Run elastic net ****"
conda run -p "${ENV_PATH}/R_env" Rscript ../_h/01.elastic-net.R

log_message "**** Job ends ****"

