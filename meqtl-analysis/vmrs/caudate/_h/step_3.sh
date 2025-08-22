#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=preprocess_genotypes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sierramannion2028@u.northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --output=genotypes.%j.log
#SBATCH --time=00:30:00

# Function to echo with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility
log_message "**** Loading modules ****"

module purge
module list

plink2 --version
## Edit with your job command
OUTDIR="protected_data"
GENOTYPES="/projects/b1213/large_projects/localQTL-benchmarking/inputs/genotypes"

log_message "**** Format genotypes ****"
mkdir -p $OUTDIR

plink2 --pfile $GENOTYPES/_m/TOPMed_LIBD \
       --keep keepPsam.txt --make-pgen \
       --out $OUTDIR/TOPMed_LIBD

log_message "**** Job ends ****"
