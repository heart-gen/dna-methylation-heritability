#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=summary_ld07_jr
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/%x.%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

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
ENV_PATH="/projects/p32505/opt/envs"

log_message "**** Job starts ****"

export ld_decay=0.7
export METHOD=joint_ridge

log_message "**** Run info ****"
echo "RUN_NAME: ${RUN_NAME}"

log_message "**** Run summary ****"
conda run -p "${ENV_PATH}/epigenomics" Rscript "../_h/02.combined-data.R"

log_message "**** Job ends ****"
