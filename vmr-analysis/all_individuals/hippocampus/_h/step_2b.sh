#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=normal       # Partition (queue) name
#SBATCH --time=06:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=16G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-24
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=res_var  # Job name
#SBATCH --output=/dev/null      # Standard output log
#SBATCH --error=/dev/null       # Standard error log

# Map chr names 
if [ "$SLURM_ARRAY_TASK_ID" -le 22 ]; then
    CHR="$SLURM_ARRAY_TASK_ID"
elif [ "$SLURM_ARRAY_TASK_ID" -eq 23 ]; then
    CHR="X"
elif [ "$SLURM_ARRAY_TASK_ID" -eq 24 ]; then
    CHR="Y"
else
    echo "Invalid SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Create log directories for each chr
LOG_DIR="logs/chr_${CHR}"
mkdir -p "$LOG_DIR"

# Redirect output and error logs to chr-specific log files
exec > >(tee -a "$LOG_DIR/res_var_${CHR}_out.log")
exec 2> >(tee -a "$LOG_DIR/res_var_${CHR}_err.log" >&2)

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
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## List current modules for reproducibility

module purge
module list

# Set path variables
ENV_PATH="/projects/p32505/opt/envs"
CPG_DIR="./cpg/chr_${CHR}"
METH_FILE="${CPG_DIR}/cpg_meth.phen"
SPLIT=5000
OUT_DIR="${CPG_DIR}/tmp_files"

mkdir -p "$OUT_DIR"

echo "Working on: Chromosome $CHR"
echo "Splitting methylation matrix into columns of $SPLIT"

# Get number of columns (based on header) and rows
NUM_COLS=$(head -1 "$METH_FILE" | tr '\t' '\n' | wc -l)
NUM_ROWS=$(wc -l < "$METH_FILE")

# A previous run may have been killed part-way through the split, which would
# leave 02b.res_var.R silently residualizing a truncated chromosome. Reuse the
# existing chunks only when the whole set is present and each has full height.
EXPECTED_CHUNKS=$(( (NUM_COLS - 2 + SPLIT - 1) / SPLIT ))
EXISTING_CHUNKS=$(ls "$OUT_DIR"/cpg_meth_*.tsv 2>/dev/null | wc -l)
SPLIT_OK=0

if [ "$EXISTING_CHUNKS" -eq "$EXPECTED_CHUNKS" ]; then
    SPLIT_OK=1
    for chunk in "$OUT_DIR"/cpg_meth_*.tsv; do
        if [ "$(wc -l < "$chunk")" -ne "$NUM_ROWS" ]; then
            echo "Truncated split file ($(wc -l < "$chunk") of $NUM_ROWS rows): $chunk"
            SPLIT_OK=0
            break
        fi
    done
fi

# Loop over column chunks
if [ "$SPLIT_OK" -eq 1 ]; then
    echo "All $EXPECTED_CHUNKS split files already exist, skipping splitting."
else
    if [ "$EXISTING_CHUNKS" -gt 0 ]; then
        echo "Discarding unusable split ($EXISTING_CHUNKS of $EXPECTED_CHUNKS files present); re-splitting."
        rm -f "$OUT_DIR"/cpg_meth_*.tsv
    fi
    for ((i=3; i<=NUM_COLS; i+=SPLIT)); do
        start=$i
        end=$((i + SPLIT - 1))
        if [ "$end" -gt "$NUM_COLS" ]; then
            end=$NUM_COLS
        fi
        echo "Extracting columns $start to $end"

        # Generate list of columns to extract (e.g., 1,2,3,...)
        cols=$(seq -s, $start $end)
        cut -f1,2,$cols "$METH_FILE" > "${OUT_DIR}/cpg_meth_${start}_${end}.tsv"
    done

    WRITTEN_CHUNKS=$(ls "$OUT_DIR"/cpg_meth_*.tsv 2>/dev/null | wc -l)
    if [ "$WRITTEN_CHUNKS" -ne "$EXPECTED_CHUNKS" ]; then
        log_message "Error: wrote $WRITTEN_CHUNKS split files, expected $EXPECTED_CHUNKS"
        exit 1
    fi
fi

## Activate conda environment
conda run -p $ENV_PATH/epigenomics Rscript ../_h/02b.res_var.R $CHR

if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
