#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=enet_5k
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/elastic_h2_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000%250
#SBATCH --time=02:00:00

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

task_id=$SLURM_ARRAY_TASK_ID
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
module list

ENV_PATH="/projects/p32505/opt/env"
export NUM_SAMPLES=5000
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

export TMPDIR="/projects/b1042/HEART-GeN-Lab/tmp"
mkdir -p $TMPDIR

log_message "**** Run elastic net ****"
conda run -p "${ENV_PATH}/r_env" Rscript ../_h/01.elastic-net.R

log_message "**** Job ends ****"
