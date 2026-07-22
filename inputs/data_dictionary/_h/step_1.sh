#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=data_dictionary
#SBATCH --output=logs/data_dictionary.%j.log

# Submit from inputs/data_dictionary/_m/:
#   mkdir -p logs && sbatch ../_h/step_1.sh

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
echo "Working directory: $(pwd)"

module purge
module list

mkdir -p logs

log_message "Building Phase 0 data dictionary tables"
python3 ../_h/00_build_data_dictionary.py

if [ $? -ne 0 ]; then
    log_message "Error: data dictionary build failed"
    exit 1
fi

log_message "**** Job ends ****"
