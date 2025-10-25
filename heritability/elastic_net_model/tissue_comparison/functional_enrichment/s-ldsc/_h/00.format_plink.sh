#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --job-name=format_plink  # Job name
#SBATCH --output=logs/format_plink/output_%j.log  # Standard output log
#SBATCH --error=logs/format_plink/error_%j.log    # Standard error log

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

echo "Processing chromosome: $chr"

# Directory containing BIM files
plink_dir="/projects/b1213/users/alexis/projects/dna-methylation-heritability/heritability/caudate/_m/plink_format/chr_${chr}"
output_bim="plink_files/output_dir/chr_${chr}.bim"

# Merge all .bim files in the directory
cat "$plink_dir"/*.bim > "$output_bim.tmp"

# Sort by the 4th column numerically
sort -k4,4n "$output_bim.tmp" > "$output_bim"

# Remove temporary file
rm "$output_bim.tmp"

echo "BIM files for chromosome $chr have been merged and sorted."

# Directory containing FAM files
output_fam="plink_files/chr_${chr}.fam"

# Merge all .bim files in the directory
cat "$plink_dir"/*.fam > "$output_fam.tmp"

# Remove temporary file
rm "$output_fam.tmp"

echo "FAM files for chromosome $chr have been merged."

log_message "**** Job ends ****"