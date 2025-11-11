#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --job-name=partition_heritability
#SBATCH --output=logs/04.output_%j.log
#SBATCH --error=logs/04.error_%j.log

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

BRAIN_REGIONS=("caudate")
HERITABILITY=("heritable_lifted")
DISEASES=("scz")

SCRIPT=../../ldsc/ldsc.py
SUMSTATS_DIR=./sumstats
BASELINE_LD_DIR=../../resource/1000G_Phase3_baselineLD_v2.2_ldscores
WEIGHTS_DIR=../../resource/1000G_Phase3_weights_hm3_no_MHC
FRQ_DIR=../../resource/1000G_Phase3_frq

for DISEASE in "${DISEASES[@]}"; do
	echo "Processing disease: ${DISEASE}"

	for REGION in "${BRAIN_REGIONS[@]}"; do

		echo "Processing region: ${REGION}"

		for STATUS in "${HERITABILITY[@]}"; do

			echo "Processing status: ${STATUS}"

			CUSTOM_LD_DIR=./custom_ldscores/${REGION}/${STATUS}
			OUT_DIR=./results/${DISEASE}/${REGION}/${STATUS}

			mkdir -p $OUT_DIR

				python $SCRIPT \
		    		--h2 $SUMSTATS_DIR/${DISEASE}.sumstats.gz \
		    		--ref-ld-chr $BASELINE_LD_DIR/baselineLD.,$CUSTOM_LD_DIR/${REGION}_${STATUS}. \
		    		--w-ld-chr $WEIGHTS_DIR/weights.hm3_noMHC. \
					--frqfile-chr $FRQ_DIR/1000G.EUR.QC. \
					--overlap-annot \
					--thin-annot \
					--print-coefficients \
		    		--out $OUT_DIR/${DISEASE}_${REGION}_${STATUS}

		done
	done
done

log_message "**** Job ends ****"
