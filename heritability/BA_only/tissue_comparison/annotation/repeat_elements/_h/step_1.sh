#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=overlap_repeat_elements
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

log_message "Overlapping VMRs with RepeatMasker repeat elements"

conda run -p /projects/p32505/opt/envs/epigenomics \
  Rscript "../_h/01.overlap_repeat_elements.R"

if [ $? -ne 0 ]; then
    log_message "Error: script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
