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

# Copy bed files
for TISSUE in "${REGIONS[@]}"; do
	cp ../../../../${TISSUE}/_m/vmr.bed ./${TISSUE}_unsorted.bed
	sort -k1,1 -k2,2n ./${TISSUE}_unsorted.bed > ./${TISSUE}.bed
done

for f in "${OVERLAP[@]}"; do
    OUTDIR="./f_${f}"
    mkdir -p "$OUTDIR"/{sets,jaccard}

	# Extract tissue specific VMRs
	bedtools intersect -a caudate.bed -b hippocampus.bed dlpfc.bed \
	-v -f $f > $OUTDIR/sets/caudate_specific.bed

	bedtools intersect -a hippocampus.bed -b caudate.bed dlpfc.bed \
	-v -f $f > $OUTDIR/sets/hippocampus_specific.bed

	bedtools intersect -a dlpfc.bed -b hippocampus.bed caudate.bed \
	-v -f $f > $OUTDIR/sets/dlpfc_specific.bed

	# Extract shared VMRs in 2 tissues
	bedtools intersect -a caudate.bed -b hippocampus.bed -f $f -wo | \
	bedtools intersect -a - -b dlpfc.bed \
	-v -f $f > $OUTDIR/sets/caudate_hippocampus_overlap.bed

	bedtools intersect -a caudate.bed -b dlpfc.bed -f $f -wo | \
	bedtools intersect -a - -b hippocampus.bed \
	-v -f $f > $OUTDIR/sets/caudate_dlpfc_overlap.bed

	bedtools intersect -a hippocampus.bed -b dlpfc.bed -f $f -wo | \
	bedtools intersect -a - -b caudate.bed \
	-v -f $f > $OUTDIR/sets/hippocampus_dlpfc_overlap.bed

	bedtools jaccard \
		-a <(bedtools intersect -a caudate.bed -b dlpfc.bed -v -f $f | sort -k1,1 -k2,2n) \
		-b <(bedtools intersect -a hippocampus.bed -b dlpfc.bed -v -f $f | sort -k1,1 -k2,2n) \
		-f $f > $OUTDIR/jaccard/caudate_hippocampus_jaccard.tsv

	bedtools jaccard \
		-a <(bedtools intersect -a caudate.bed -b hippocampus.bed -v -f $f | sort -k1,1 -k2,2n) \
		-b <(bedtools intersect -a dlpfc.bed -b hippocampus.bed -v -f $f | sort -k1,1 -k2,2n) \
		-f $f > $OUTDIR/jaccard/caudate_dlpfc_jaccard.tsv

	bedtools jaccard \
		-a <(bedtools intersect -a hippocampus.bed -b caudate.bed -v -f $f | sort -k1,1 -k2,2n) \
		-b <(bedtools intersect -a dlpfc.bed -b caudate.bed -v -f $f | sort -k1,1 -k2,2n) \
		-f $f > $OUTDIR/jaccard/hippocampus_dlpfc_jaccard.tsv

	# Extract shared VMRs in all tissues
	bedtools intersect -a caudate.bed -b dlpfc.bed -f $f -wo | \
	bedtools intersect -a - -b hippocampus.bed \
	-f $f -wo > $OUTDIR/sets/3tissues_overlap.bed

	bedtools jaccard \
		-a <(bedtools intersect -a caudate.bed -b hippocampus.bed -f $f | sort -k1,1 -k2,2n) \
		-b dlpfc.bed \
		-f $f > $OUTDIR/jaccard/3tissues_jaccard.tsv

    if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
    fi
done

for F in "${OVERLAP[@]}"; do
    OUTDIR="./F_${F}"
    mkdir -p "$OUTDIR"/{sets,jaccard,percent_overlap}    

	# Extract tissue specific VMRs
	bedtools intersect -a caudate.bed -b hippocampus.bed dlpfc.bed \
	-v -F $F > $OUTDIR/sets/caudate_specific.bed

	bedtools intersect -a hippocampus.bed -b caudate.bed dlpfc.bed \
	-v -F $F > $OUTDIR/sets/hippocampus_specific.bed

	bedtools intersect -a dlpfc.bed -b hippocampus.bed caudate.bed \
	-v -F $F > $OUTDIR/sets/dlpfc_specific.bed

	# Extract shared VMRs in 2 tissues
	bedtools intersect -a caudate.bed -b hippocampus.bed -F $F -wo | \
	bedtools intersect -a - -b dlpfc.bed \
	-v -F $F > $OUTDIR/sets/caudate_hippocampus_overlap.bed

	bedtools intersect -a caudate.bed -b dlpfc.bed -F $F -wo | \
	bedtools intersect -a - -b hippocampus.bed \
	-v -F $F > $OUTDIR/sets/caudate_dlpfc_overlap.bed

	bedtools intersect -a hippocampus.bed -b dlpfc.bed -F $F -wo | \
	bedtools intersect -a - -b caudate.bed \
	-v -F $F > $OUTDIR/sets/hippocampus_dlpfc_overlap.bed

	bedtools jaccard \
		-a <(bedtools intersect -a caudate.bed -b dlpfc.bed -v -F $F | sort -k1,1 -k2,2n) \
		-b <(bedtools intersect -a hippocampus.bed -b dlpfc.bed -v -F $F | sort -k1,1 -k2,2n) \
		-F $F > $OUTDIR/jaccard/caudate_hippocampus_jaccard.tsv

	bedtools jaccard \
		-a <(bedtools intersect -a caudate.bed -b hippocampus.bed -v -F $F | sort -k1,1 -k2,2n) \
		-b <(bedtools intersect -a dlpfc.bed -b hippocampus.bed -v -F $F | sort -k1,1 -k2,2n) \
		-F $F > $OUTDIR/jaccard/caudate_dlpfc_jaccard.tsv

	bedtools jaccard \
		-a <(bedtools intersect -a hippocampus.bed -b caudate.bed -v -F $F | sort -k1,1 -k2,2n) \
		-b <(bedtools intersect -a dlpfc.bed -b caudate.bed -v -F $F | sort -k1,1 -k2,2n) \
		-F $F > $OUTDIR/jaccard/hippocampus_dlpfc_jaccard.tsv

	# Extract shared VMRs in all tissues
	bedtools intersect -a caudate.bed -b dlpfc.bed -F $F -wo | \
	bedtools intersect -a - -b hippocampus.bed \
	-F $F -wo > $OUTDIR/sets/3tissues_overlap.bed

	bedtools jaccard \
		-a <(bedtools intersect -a caudate.bed -b hippocampus.bed -F $F | sort -k1,1 -k2,2n) \
		-b dlpfc.bed \
		-F $F > $OUTDIR/jaccard/3tissues_jaccard.tsv

    if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
    fi
done

#conda deactivate

log_message "**** Job ends ****"
