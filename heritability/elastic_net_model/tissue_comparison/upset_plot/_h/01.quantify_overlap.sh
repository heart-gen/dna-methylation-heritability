#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=00:10:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=8G                # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=overlap # Job name
#SBATCH --output=logs/overlap.%j.log # Standard output log

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

## List current modules for reproducibility

module purge
module load bedtools/2.31.1
module list

# Set path variables
log_message "Calculating brain region overlap"
#conda activate /projects/p32505/opt/envs/genomics
#export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

REGIONS=(caudate hippocampus dlpfc)
OVERLAP=(0.25 0.5 0.75)

## Function
extract_tissue_specific(){
	local tissue1=$1 tissue2=$2 tissue3=$3 overlap=$4 out_file=$5
	bedtools intersect -a "$tissue1" -b "$tissue2" "$tissue3" \
	"$FLAG" "$overlap" -v > "$out_file"
}

extract_pairwise_overlap(){
	local tissue1=$1 tissue2=$2 overlap=$3 out_file=$4
	bedtools intersect -a "$tissue1" -b "$tissue2" \
	"$FLAG" "$overlap" -wo > "$out_file"
}

extract_3tissue(){
	local tissue1=$1 tissue2=$2 tissue3=$3 overlap=$4 out_bed=$5 out_jaccard=$6
	local tmp_bed="${out_bed}.tmp"

	bedtools intersect -a "$tissue1" -b "$tissue2" "$FLAG" "$overlap" -u | \
	bedtools intersect -a - -b "$tissue3" "$FLAG" "$overlap" -u > "$tmp_bed"

	bedtools multiinter \
		-header -names "$tissue1" "$tissue2" "$tissue3" \
		-i "$tmp_bed" "$tissue2" "$tissue3" > "$out_bed"

	bedtools jaccard \
		-a <(bedtools intersect -a "$tissue1" -b "$tissue2" \
		"$FLAG" "$overlap" | sort -k1,1 -k2,2n) \
		-b "$tissue3" "$FLAG" "$overlap" > "$out_jaccard"
}

jaccard(){
	local tissue1=$1 tissue2=$2 overlap=$3 out_file=$4
	bedtools jaccard -a "$tissue1" -b "$tissue2" \
	"$FLAG" "$overlap" > "$out_file"
}

## Main
# Copy bed files
for TISSUE in "${REGIONS[@]}"; do
	BED="./${TISSUE}.bed"
	if [ ! -f "$BED" ]; then
		log_message "Sorting bed file for ${TISSUE}."
		cp ../../../../${TISSUE}/_m/vmr.bed ./${TISSUE}_unsorted.bed
		sort -k1,1 -k2,2n ./${TISSUE}_unsorted.bed > ./${TISSUE}.bed
	else
		log_message "Bed file already exists. Skipping."
	fi
done

for f in "${OVERLAP[@]}"; do
	for FLAG in "-f" "-F"; do 
		OUTDIR="./${FLAG//-/}_${f}"
		mkdir -p "$OUTDIR"/{sets,jaccard}

		# Extract tissue specific VMRs
		extract_tissue_specific caudate.bed hippocampus.bed dlpfc.bed \
		"$f" "$OUTDIR/sets/caudate_specific.bed"
		extract_tissue_specific hippocampus.bed caudate.bed dlpfc.bed \
		"$f" "$OUTDIR/sets/hippocampus_specific.bed"
		extract_tissue_specific dlpfc.bed caudate.bed hippocampus.bed \
		"$f" "$OUTDIR/sets/dlpfc_specific.bed"

		# Extract shared VMRs in 2 tissues
		extract_pairwise_overlap caudate.bed hippocampus.bed \
		"$f" "$OUTDIR/sets/caudate_hippocampus_overlap.bed"
		extract_pairwise_overlap caudate.bed dlpfc.bed \
		"$f" "$OUTDIR/sets/caudate_dlpfc_overlap.bed"
		extract_pairwise_overlap hippocampus.bed dlpfc.bed \
		"$f" "$OUTDIR/sets/hippocampus_dlpfc_overlap.bed"

		# Calculate pairwise jaccard
		jaccard caudate.bed hippocampus.bed \
		"$f" "$OUTDIR/jaccard/caudate_hippocampus_jaccard.tsv"
		jaccard caudate.bed dlpfc.bed \
		"$f" "$OUTDIR/jaccard/caudate_dlpfc_jaccard.tsv"
		jaccard hippocampus.bed dlpfc.bed \
		"$f" "$OUTDIR/jaccard/hippocampus_dlpfc_jaccard.tsv"

		# Extract shared VMRs in all tissues
		extract_3tissue caudate.bed hippocampus.bed dlpfc.bed \
		"$f" "$OUTDIR/sets/3tissues_overlap.bed" \
		"$OUTDIR/jaccard/3tissues_jaccard.tsv"

		if [ $? -ne 0 ]; then
		log_message "Error: Conda or script execution failed"
		exit 1
		fi
	done
done

#conda deactivate

log_message "**** Job ends ****"
