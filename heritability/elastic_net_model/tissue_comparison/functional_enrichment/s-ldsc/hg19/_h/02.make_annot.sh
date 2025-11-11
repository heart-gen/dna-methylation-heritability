#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name=make_annot
#SBATCH --output=logs/02.output_make_annot.log
#SBATCH --error=logs/02.error_make_annot_.log

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

module load bedtools/2.30.0

BRAIN_REGIONS=("caudate")
HERITABILITY=("heritable_lifted")

SCRIPT_DIR=../../ldsc
BIM_DIR=../../resource/1000G_Phase3_plinkfiles

for REGION in "${BRAIN_REGIONS[@]}"; do

  echo "Processing region: $REGION"

  for STATUS in "${HERITABILITY[@]}"; do
    echo "Processing status: $STATUS"

    BED_FILE=./vmr/${REGION}/${STATUS}.bed
    OUT_DIR=./custom_ldscores/${REGION}/${STATUS}
    mkdir -p $OUT_DIR

    for CHR in {1..22}; do
      echo "Processing chromosome $CHR..."

      BIM_FILE=${BIM_DIR}/1000G.EUR.QC.${CHR}.bim

      python $SCRIPT_DIR/make_annot.py \
        --bed-file ${BED_FILE} \
        --bimfile ${BIM_FILE} \
        --annot-file $OUT_DIR/${REGION}_${STATUS}.${CHR}.annot.gz \
        --windowsize 500000 \
        
		done
  done
done

log_message "**** Job ends ****"