#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_data_audit
#SBATCH --output=logs/data_audit.%j.log

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

log_message "Auditing meQTL validation inputs"
python3 ../_h/00_data_audit.py

if [ $? -ne 0 ]; then
    log_message "Error: data audit failed"
    exit 1
fi

log_message "**** Job ends ****"
