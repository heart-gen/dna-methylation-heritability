#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name=munge_sumstats
#SBATCH --output=logs/01.output_%j.log
#SBATCH --error=logs/01.error_%j.log

# Log function
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

SCRIPT=../../ldsc/munge_sumstats.py
GWAS=../../gwas/ad_gwas.txt.gz
OUT_DIR=./sumstats
mkdir -p $OUT_DIR

python $SCRIPT \
    --sumstats $GWAS \
    --out $OUT_DIR/ad \
    --a1 effect_allele \
    --a2 other_allele \
    --signed-sumstats beta,0 \
    --p p_value \
    --snp rsID \
    --N 1126563

log_message "**** Job ends ****"