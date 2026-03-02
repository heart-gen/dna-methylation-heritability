#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=20gb
#SBATCH --job-name=run_susie
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/run_susie.%j.log

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

conda activate /projects/p32505/opt/envs/epigenomics

OUT_DIR=./results/caudate/chr_1
mkdir -p $OUT_DIR

/projects/p32505/opt/envs/epigenomics/bin/Rscript ../_h/susie.R \
  ./gwas/pd/caudate/chr_1/TOPMed_LIBD.AA.903969_904084_gwas_subset.tsv \
  ./ld_matrices/caudate/chr_1/TOPMed_LIBD.AA.903969_904084_final.ld \
  $OUT_DIR/TOPMed_LIBD.AA.903969_904084_results \

conda deactivate

log_message "**** Job ends ****"