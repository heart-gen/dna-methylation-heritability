#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=preprocess_genotypes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --output=logs/genotypes.%j.log
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
module load plink/2.0-alpha-3.3
module list

plink2 --version
## Edit with your job command
OUTDIR="./protected_data"
GENOTYPES="../../../../inputs/genotypes"
SAMPLES="../../../../heritability/caudate/_m/samples.txt"

log_message "**** Format genotypes ****"
mkdir -p $OUTDIR

plink2 --pfile $GENOTYPES/TOPMed_LIBD \
       --keep $SAMPLES --make-pgen \
       --no-parents \
       --no-sex \
       --no-pheno \
       --out $OUTDIR/TOPMed_LIBD

log_message "**** Job ends ****"
