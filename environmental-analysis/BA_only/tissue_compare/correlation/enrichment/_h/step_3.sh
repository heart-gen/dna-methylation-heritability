#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=pub_enrichment
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/pub_enrichment_%j.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

mkdir -p logs publication_figures

log_message "**** Job starts ****"

echo "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME:-N/A}"
echo "Hostname: ${HOSTNAME}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module list

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh

log_message "**** Plotting publication enrichment figures ****"
conda activate /projects/p32505/opt/envs/epigenomics
Rscript "../_h/03.plot_publication_figures.R"
conda deactivate

log_message "**** Job ends ****"
