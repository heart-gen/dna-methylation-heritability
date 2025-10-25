#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=1G                # Memory limit
#SBATCH --job-name=02.multi_core  # Job name
#SBATCH --output=output_%j.log  # Standard output log
#SBATCH --error=error_%j.log    # Standard error log

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

output_dir=$(dirname "$output_prefix")
mkdir -p "$output_dir"

# Sort bed file by chromosome and range positions
input_file="/projects/b1213/users/alexis/projects/dna-methylation-heritability/heritability/caudate/_m/vmr.bed"
output_prefix="vmr"

sort -k1,1 -k2,2n "$input_file" | awk '{print > "'"${output_prefix}"'_chr"$1".bed"}'

log_message "**** Job ends ****"