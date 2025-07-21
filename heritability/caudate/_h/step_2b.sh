#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=10G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-22
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=res_var  # Job name
#SBATCH --output=/dev/null      # Standard output log
#SBATCH --error=/dev/null       # Standard error log

# Create log directories for each chr
LOG_DIR="logs/chr_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$LOG_DIR"

# Redirect output and error logs to chr-specific log files
exec > >(tee -a "$LOG_DIR/res_var_${SLURM_ARRAY_TASK_ID}_out.log")
exec 2> >(tee -a "$LOG_DIR/res_var_${SLURM_ARRAY_TASK_ID}_err.log" >&2)

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
ENV_PATH="/projects/p32505/opt/env"
CPG_DIR="./cpg/chr_${SLURM_ARRAY_TASK_ID}"
METH_FILE="${CPG_DIR}/cpg_meth.phen"
SPLIT=5000
OUT_DIR="${CPG_DIR}/tmp_files"

mkdir -p "$OUT_DIR"

echo "Working on: Chromosome $SLURM_ARRAY_TASK_ID"
echo "Splitting methylation matrix into columns of $SPLIT"

# Get number of columns (based on header)
NUM_COLS=$(head -1 "$METH_FILE" | tr '\t' '\n' | wc -l)

# Loop over column chunks
if ls "$OUT_DIR"/cpg_meth_*.tsv 1> /dev/null 2>&1; then
    echo "Split files already exist, skipping splitting."
else
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
fi

## Activate conda environment
conda run -p $ENV_PATH/r_env Rscript ../_h/02b.res_var.R $SLURM_ARRAY_TASK_ID

if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

log_message "**** Job ends ****"
