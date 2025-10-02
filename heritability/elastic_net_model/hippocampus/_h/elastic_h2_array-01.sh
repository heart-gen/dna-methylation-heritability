#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=elastic_h2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --output=logs/elastic_h2_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --array=1-9576%300
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

task_id=${SLURM_ARRAY_TASK_ID}
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
module list

ENV_PATH="/projects/p32505/opt/envs"
export region="hippocampus"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

log_message "**** Run elastic net ****"
conda run -p "${ENV_PATH}/epigenomics" Rscript ../_h/01.elastic-net.R

log_message "**** Job ends ****"

