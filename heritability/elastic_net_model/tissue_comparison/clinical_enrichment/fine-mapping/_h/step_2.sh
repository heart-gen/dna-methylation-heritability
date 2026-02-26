#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=02:30:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --job-name=subset_gwas  # Job name
#SBATCH --output=logs/step_2/output_%j.log  # Standard output log
#SBATCH --error=logs/step_2/error_%j.log    # Standard error log
#SBATCH --array=1-22

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

module load plink/1.9

REF_DIR=/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/caudate/_m/plink_format/chr_$SLURM_ARRAY_TASK_ID
GWAS_FILE="/projects/b1213/resources/gwas/mdd/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.gz"

# Loop through each BED file in plink directory
for file in "$REF_DIR"/TOPMed_LIBD.AA.*.bed; do
	filename=$(basename "$file")

	# Extract START and END positions from filename
	if [[ $filename =~ ([0-9]+)_([0-9]+)\. ]]; then
		START=${BASH_REMATCH[1]}
		END=${BASH_REMATCH[2]}
		echo "File: $filename → START=$START, END=$END"
	else
		echo "File $filename does not match the expected pattern."
	fi

  # Create output directory for LD matrices
	OUT_DIR=gwas/mdd/caudate/chr_$SLURM_ARRAY_TASK_ID
	mkdir -p $OUT_DIR

  # Subset GWAS file based on SNPs in the BIM file
  zcat "$GWAS_FILE" | \
  awk -v chr="$SLURM_ARRAY_TASK_ID" -v start="$START - 500000" -v end="$END + 500000" '
  NR==1 || ($2==chr && $3>=start && $3<=end)
  ' > "$OUT_DIR/TOPMed_LIBD.AA.${START}_${END}.tsv"
done

log_message "**** Job ends ****"