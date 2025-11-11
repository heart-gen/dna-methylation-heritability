#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --array=1-22
#SBATCH --job-name=compute_ldscores
#SBATCH --output=logs/03.output_%a.log
#SBATCH --error=logs/03.error_%a.log

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
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

BRAIN_REGIONS=("caudate")
HERITABILITY=("heritable_lifted")

SCRIPT=../../ldsc/ldsc.py
BIM_DIR=../../resource/1000G_Phase3_plinkfiles
HAPMAP3_SNPS=../../resource/hm3_no_MHC.list.txt

for REGION in "${BRAIN_REGIONS[@]}"; do
	for STATUS in "${HERITABILITY[@]}"; do

    	OUT_DIR=./custom_ldscores/${REGION}/${STATUS}

    	python $SCRIPT \
			--l2 \
			--bfile $BIM_DIR/1000G.EUR.QC.${SLURM_ARRAY_TASK_ID} \
			--ld-wind-cm 1 \
			--annot $OUT_DIR/${REGION}_${STATUS}.${SLURM_ARRAY_TASK_ID}.annot.gz \
			--thin-annot \
			--out $OUT_DIR/${REGION}_${STATUS}.${SLURM_ARRAY_TASK_ID} \
			--print-snps $HAPMAP3_SNPS

	done
done

log_message "**** Job ends ****"