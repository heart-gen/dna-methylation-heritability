#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=clean_ld07
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=%x_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_ROOT="${PROJECT_DIR}/_m"
RUN_NAME="ld_0.7_sim_200_indiv"
RUN_DIR="${OUTPUT_ROOT}/${RUN_NAME}"
LOG_DIR="${RUN_DIR}/logs"

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG_DIR}/clean_data_${SLURM_JOBID:-manual}.log")
exec 2>&1

log_message "**** Job starts ****"

num_samples=200
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

log_message "**** Cleaning directory ****"
gzip -9v "simulation_${num_samples}_h2_elastic-net.tsv"
gzip -9v "simulation_${num_samples}_betas_elastic-net.tsv"

rm -r betas/ summary/ h2/

log_message "**** Job ends ****"
