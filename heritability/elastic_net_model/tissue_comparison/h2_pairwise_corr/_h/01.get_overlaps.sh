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
log_message "Writing pairwise brain region overlap"
#conda activate /projects/p32505/opt/envs/genomics
#export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

OVERLAP=(0.25 0.5 0.75)

for F in "${OVERLAP[@]}"; do
    OUTDIR="./F_${F}"
    mkdir -p "$OUTDIR"

	# Counts
	bedtools intersect -a ../../../../caudate/_m/vmr.bed \
	-b ../../../../hippocampus/_m/vmr.bed -F $F -wa -wb | \
	bedtools intersect -a - -b ../../../../dlpfc/_m/vmr.bed \
	-v -F $F > $OUTDIR/caudate_hippocampus_overlap_${F}.bed

	bedtools intersect -a ../../../../caudate/_m/vmr.bed \
	-b ../../../../dlpfc/_m/vmr.bed -F $F -wa -wb | \
	bedtools intersect -a - -b ../../../../hippocampus/_m/vmr.bed \
	-v -F $F > $OUTDIR/caudate_dlpfc_overlap_${F}.bed

	bedtools intersect -a ../../../../hippocampus/_m/vmr.bed \
	-b ../../../../dlpfc/_m/vmr.bed -F $F -wa -wb | \
	bedtools intersect -a - -b ../../../../caudate/_m/vmr.bed \
	-v -F $F > $OUTDIR/hippocampus_dlpfc_overlap_${F}.bed

    if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
    fi
done

#conda deactivate

log_message "**** Job ends ****"
