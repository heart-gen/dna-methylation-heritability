#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --mem=10G
#SBATCH --job-name=elastic_h2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --output=logs/elastic_h2_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --array=1-11575%250
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

task_id=${SLURM_ARRAY_TASK_ID}
export task_id

echo "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "OFFSET: ${OFFSET}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Computed task_id: ${task_id}"

module purge
module load plink/2.0-alpha-3.3
module list

ENV_PATH="/projects/p32505/opt/envs"
export region="caudate"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

log_message "**** Run elastic net ****"
conda run -p "${ENV_PATH}/epigenomics" Rscript ../_h/01.elastic-net.R

log_message "**** Job ends ****"
