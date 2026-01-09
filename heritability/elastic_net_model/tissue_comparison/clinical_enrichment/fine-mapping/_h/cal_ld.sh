#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=run
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=logs/cal_ld.%A-%a.log
#SBATCH --array=432
#SBATCH --time=01:00:00

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

REGION_LIST="/projects/p32505/users/elisa/dna-methylation-heritability/heritability/caudate/_m/vmr.bed"

# Get the current region name from the region list
REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $REGION_LIST)
CHR=$(echo "$REGION" | awk '{print $1}')
START=$(echo "$REGION" | awk '{print $2}')
START_POSITION=$((START - 500000))
END=$(echo "$REGION" | awk '{print $3}')
END_POSITION=$((END + 500000))

python ../_h/SuSiEx_LD.py \
	--ref_file=AFR \
	--ld_file=AFR \
	--chr=$CHR \
	--bp=$START_POSITION,$END_POSITION \
	--plink=/projects/b1213/users/alexis/projects/dna-methylation-heritability/heritability/caudate/_m/plink_format/chr_$CHR/TOPMed_LIBD.AA.${START}_${END} \
	--maf=0.005

log_message "**** Job ends ****"