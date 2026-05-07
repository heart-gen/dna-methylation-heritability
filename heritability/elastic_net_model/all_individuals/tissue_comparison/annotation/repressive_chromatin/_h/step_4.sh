#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=plot_region_specificity
#SBATCH --output=logs/%x.%j.log

set -euo pipefail


log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module list

mkdir -p logs

log_message "Plotting region-specificity figures for repressive heritable VMRs"

conda run -p /projects/p32505/opt/envs/epigenomics \
  Rscript "../_h/04.plot_region_specificity.R"

if [ $? -ne 0 ]; then
    log_message "Error: script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
