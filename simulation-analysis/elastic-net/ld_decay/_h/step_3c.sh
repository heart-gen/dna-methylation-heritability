#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=clean_ld07
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
RUN_NAME="ld_0.7_sim_200_indiv"

log_message "**** Job starts ****"

ld_decay=0.7

log_message "**** Run info ****"
echo "RUN_NAME: ${RUN_NAME}"

log_message "**** Cleaning directory ****"
for method in boosting_hybrid joint_ridge; do
    gzip -9v "simulation_200_${ld_decay}_${method}_h2_elastic-net.tsv"
    gzip -9v "simulation_200_${ld_decay}_${method}_betas_elastic-net.tsv"
    rm -r betas_${method}/ summary_${method}/ h2_${method}/
done

log_message "**** Job ends ****"
