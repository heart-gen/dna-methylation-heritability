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
module list

# Set path variables
log_message "Calculating brain region overlap"
conda activate /projects/p32505/opt/envs/genomics
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

REGIONS=(caudate hippocampus dlpfc)
OVERLAP=(0.25 0.5 0.75)

# Extract tissue specific
for f in "${OVERLAP[@]}"; do
    OUTDIR="./f_${f}"
    mkdir -p "$OUTDIR"/{sets,jaccard,percent_overlap}

	# Counts
	bedtools intersect -a caudate.bed -b hippocampus.bed dlpfc.bed \
	-v -f $f > $OUTDIR/sets/caudate_specific_${f}.bed

	bedtools intersect -a hippocampus.bed -b caudate.bed dlpfc.bed \
	-v -f $f > $OUTDIR/sets/hippocampus_specific_${f}.bed

	bedtools intersect -a dlpfc.bed -b hippocampus.bed dlpfc.bed \
	-v -f $f > $OUTDIR/sets/dlpfc_specific_${f}.bed

    if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
    fi
done

for F in "${OVERLAP[@]}"; do
    OUTDIR="./f_${F}"
    mkdir -p "$OUTDIR"/{sets,jaccard,percent_overlap}    

	# Counts
	bedtools intersect -a caudate.bed -b hippocampus.bed dlpfc.bed \
	-v -F $F > $OUTDIR/sets/caudate_specific_${F}.bed

	bedtools intersect -a hippocampus.bed -b caudate.bed dlpfc.bed \
	-v -F $F > $OUTDIR/sets/hippocampus_specific_${F}.bed

	bedtools intersect -a dlpfc.bed -b hippocampus.bed dlpfc.bed \
	-v -F $F > $OUTDIR/sets/dlpfc_specific_${F}.bed

    if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
    fi
done

# Extract shared VMRs in 2 tissues
	bedtools intersect -a caudate.bed -b hippocampus.bed \
	-u -F $F > $OUTDIR/sets/caudate_hippocampus_overlap_${F}.bed

	bedtools jaccard -a caudate.bed -b hippocampus.bed \
	-F $F > $OUTDIR/jaccard/caudate_hippocampus_jaccard_${F}.bed

# Extract shared VMRs in all tissues
	bedtools intersect -a dlpfc.bed -b hippocampus.bed dlpfc.bed \
	-u -F $F > $OUTDIR/sets/3tissues_overlap_${F}.bed

	bedtools jaccard -a caudate.bed -b hippocampus.bed dlpfc.bed \
	-F $F > $OUTDIR/jaccard/3tissues_jaccard_${F}.bed

conda deactivate

log_message "**** Job ends ****"
