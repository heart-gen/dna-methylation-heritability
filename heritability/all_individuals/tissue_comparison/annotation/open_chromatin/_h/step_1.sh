#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=vmr_atac_overlap
#SBATCH --output=logs/vmr_atac_overlap.%j.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module list

ENV_PATH="/projects/p32505/opt/envs"

mkdir -p logs

log_message "Overlapping intergenic VMRs with brain cell-type open chromatin peaks"

conda run -p $ENV_PATH/epigenomics Rscript ../_h/01.overlap_open_chromatin.R

if [ $? -ne 0 ]; then
    log_message "Error: script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
