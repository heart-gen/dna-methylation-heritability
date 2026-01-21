#!/bin/bash
#SBATCH --account=b1213        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=1G                # Memory limit
#SBATCH --job-name=cal_ld  # Job name
#SBATCH --output=logs/output_%j.log  # Standard output log
#SBATCH --error=logs/error_%j.log    # Standard error log
#SBATCH --array=1-22         # Array job range (if needed)

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

mkdir -p cal_ld
cd cal_ld
mkdir -p chr_${SLURM_ARRAY_TASK_ID}
cd ..

REF_DIR=/projects/b1213/users/alexis/projects/dna-methylation-heritability/heritability/caudate/_m/plink_format/chr_${SLURM_ARRAY_TASK_ID}

for file in "$REF_DIR"/TOPMed_LIBD.AA.*.bed; do
	filename=$(basename "$file")

	if [[ $filename =~ ([0-9]+)_([0-9]+)\. ]]; then
		START=${BASH_REMATCH[1]}
		END=${BASH_REMATCH[2]}
		echo "File: $filename → START=$START, END=$END"
	else
		echo "File $filename does not match the expected pattern."
	fi

	mkdir -p cal_ld/chr_${SLURM_ARRAY_TASK_ID}/${START}_${END}

	START_POSITION=$((START - 500000))
	END_POSITION=$((END + 500000))

	python ../SuSiEx/utilities/SuSiEx_LD.py \
		--ref_file=$REF_DIR/TOPMed_LIBD.AA.${START}_${END} \
		--ld_file=cal_ld/chr_${SLURM_ARRAY_TASK_ID}/${START}_${END}/TOPMed_LIBD.AA.${START}_${END} \
		--chr=${SLURM_ARRAY_TASK_ID} \
		--bp=$START_POSITION,$END_POSITION \
		--plink=../SuSiEx/utilities/plink \
		--maf=0.005

done

log_message "**** Job ends ****"