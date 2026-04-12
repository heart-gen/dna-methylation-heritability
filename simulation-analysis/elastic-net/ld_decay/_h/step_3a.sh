#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=clean_ld05
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/%x.%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

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

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG_DIR}/clean_data_${SLURM_JOBID:-manual}.log")
exec 2>&1

log_message "**** Job starts ****"

ld_decay=0.5
mkdir -p "${RUN_DIR}"

log_message "**** Run info ****"
echo "Computed task_id: ${task_id}"
echo "RUN_NAME: ${RUN_NAME}"
echo "RUN_DIR: ${RUN_DIR}"

log_message "**** Cleaning directory ****"
gzip -9v "simulation_200_${ld_decay}_h2_elastic-net.tsv"
gzip -9v "simulation_200_${ld_decay}_betas_elastic-net.tsv"

rm -r betas/ summary/ h2/

log_message "**** Job ends ****"
