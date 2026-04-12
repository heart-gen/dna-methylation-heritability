#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=summary_ld05
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=%x_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_ROOT="${PROJECT_DIR}/_m"
RUN_NAME="ld_0.5_sim_200_indiv"
RUN_DIR="${OUTPUT_ROOT}/${RUN_NAME}"
LOG_DIR="${RUN_DIR}/logs"
ENV_PATH="/projects/p32505/opt/envs"

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG_DIR}/summary_data_${SLURM_JOBID:-manual}.log")
exec 2>&1

log_message "**** Job starts ****"

export num_samples=200
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

echo "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "num_samples: ${num_samples}"
echo "RUN_DIR: ${RUN_DIR}"

module purge
module list

log_message "**** Run summary ****"
conda run -p "${ENV_PATH}/epigenomics" Rscript "${SCRIPT_DIR}/02.combined-data.R"

log_message "**** Job ends ****"
