#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=load_caudate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/load_caudate.%j.log
#SBATCH --time=18:00:00
#SBATCH --mem=100G

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"

log_message "**** Loading modules ****"
module purge
module list

log_message "**** Loading conda environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/epigenomics

log_message "**** Running caudate data loading ****"
Rscript ../_h/01.load_bsseq_caudate.R

status=$?
if [ $status -ne 0 ]; then
    log_message "Error: R execution failed (status=$status)"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
