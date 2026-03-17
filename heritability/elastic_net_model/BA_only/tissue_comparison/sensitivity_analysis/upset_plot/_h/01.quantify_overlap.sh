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
BED_DIR="./bed_files"
for f in "${OVERLAP[@]}"; do
	for FLAG in "-f" "-F"; do 
		for H2_CAT in "all" "heritable" "non-heritable" "low_prediction"; do
			OUTDIR="./${FLAG//-/}_${f}"
			mkdir -p "$OUTDIR"/{sets,jaccard}

			# Extract tissue specific VMRs
			extract_tissue_specific $BED_DIR/caudate_${H2_CAT}.bed $BED_DIR/hippocampus_${H2_CAT}.bed $BED_DIR/dlpfc_${H2_CAT}.bed \
			"$f" "$OUTDIR/sets/caudate_specific_${H2_CAT}.bed"
			extract_tissue_specific $BED_DIR/hippocampus_${H2_CAT}.bed $BED_DIR/caudate_${H2_CAT}.bed $BED_DIR/dlpfc_${H2_CAT}.bed \
			"$f" "$OUTDIR/sets/hippocampus_specific_${H2_CAT}.bed"
			extract_tissue_specific $BED_DIR/dlpfc_${H2_CAT}.bed $BED_DIR/caudate_${H2_CAT}.bed $BED_DIR/hippocampus_${H2_CAT}.bed \
			"$f" "$OUTDIR/sets/dlpfc_specific_${H2_CAT}.bed"

			# Extract shared VMRs in 2 tissues
			extract_pairwise_overlap $BED_DIR/caudate_${H2_CAT}.bed $BED_DIR/hippocampus_${H2_CAT}.bed \
			"$f" "$OUTDIR/sets/caudate_hippocampus_overlap_${H2_CAT}.bed"
			extract_pairwise_overlap $BED_DIR/caudate_${H2_CAT}.bed $BED_DIR/dlpfc_${H2_CAT}.bed \
			"$f" "$OUTDIR/sets/caudate_dlpfc_overlap_${H2_CAT}.bed"
			extract_pairwise_overlap $BED_DIR/hippocampus_${H2_CAT}.bed $BED_DIR/dlpfc_${H2_CAT}.bed \
			"$f" "$OUTDIR/sets/hippocampus_dlpfc_overlap_${H2_CAT}.bed"

			# Calculate pairwise jaccard
			jaccard $BED_DIR/caudate_${H2_CAT}.bed $BED_DIR/hippocampus_${H2_CAT}.bed \
			"$f" "$OUTDIR/jaccard/caudate_hippocampus_jaccard_${H2_CAT}.tsv"
			jaccard $BED_DIR/caudate_${H2_CAT}.bed $BED_DIR/dlpfc_${H2_CAT}.bed \
			"$f" "$OUTDIR/jaccard/caudate_dlpfc_jaccard_${H2_CAT}.tsv"
			jaccard $BED_DIR/hippocampus_${H2_CAT}.bed $BED_DIR/dlpfc_${H2_CAT}.bed \
			"$f" "$OUTDIR/jaccard/hippocampus_dlpfc_jaccard_${H2_CAT}.tsv"

			# Extract shared VMRs in all tissues
			extract_3tissue $BED_DIR/caudate_${H2_CAT}.bed $BED_DIR/hippocampus_${H2_CAT}.bed $BED_DIR/dlpfc_${H2_CAT}.bed \
			"$f" "$OUTDIR/sets/3tissues_overlap_${H2_CAT}.bed" \
			"$OUTDIR/jaccard/3tissues_jaccard_${H2_CAT}.tsv"

			if [ $? -ne 0 ]; then
			log_message "Error: Conda or script execution failed"
			exit 1
			fi
		done
	done
done

#conda deactivate

log_message "**** Job ends ****"
