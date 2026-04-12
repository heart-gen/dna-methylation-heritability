#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=enet_ld05_jr
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/%x.%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000%200
#SBATCH --time=02:00:00

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

echo "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME:-N/A}"
echo "Hostname: ${HOSTNAME}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module list

SCRIPT_DIR="../_h"
OUTPUT_ROOT="../_m"
RUN_NAME="ld_0.5_sim_200_indiv"
RUN_DIR="${OUTPUT_ROOT}/${RUN_NAME}"
LOG_DIR="${RUN_DIR}/logs"

ENV_PATH="/projects/p32505/opt/envs"

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG_DIR}/elastic_h2.${SLURM_JOBID:-manual}_${SLURM_ARRAY_TASK_ID:-na}.log")
exec 2>&1

log_message "**** Job starts ****"

task_id=${SLURM_ARRAY_TASK_ID}
export task_id
export RUN_NAME
export METHOD=joint_ridge
export NUM_SAMPLES=200
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

export TMPDIR="/projects/b1042/HEART-GeN-Lab/tmp"
mkdir -p "${TMPDIR}"

mkdir -p "${RUN_DIR}"

log_message "**** Run info ****"
echo "Computed task_id: ${task_id}"
echo "RUN_NAME: ${RUN_NAME}"
echo "RUN_DIR: ${RUN_DIR}"

log_message "**** Run elastic net ****"
conda run -p "${ENV_PATH}/epigenomics" Rscript "${SCRIPT_DIR}/01.elastic-net.R"

log_message "**** Job ends ****"
