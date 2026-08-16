#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=make_figures
#SBATCH --mail-type=FAIL ## If you want to have it email you for any reason
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu ## replace with your email
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --output=logs/make_figures.%j.log
#SBATCH --time=01:00:00

# Function to echo with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility
log_message "**** Loading modules ****"

module purge
module list

# Set path variables
log_message "**** Loading conda environment ****"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

python ../_h/fdr_correction.py
rm combined_Coeff_zscore.tsv

python ../_h/make_heatmap.py

conda deactivate

log_message "**** Job ends ****"