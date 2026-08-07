#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=scz_define_loci
#SBATCH --output=logs/scz_define_loci.%j.log

# Submit from meqtl-validation/08_schizophrenia_risk_application/_m/:
#   mkdir -p logs && sbatch --export=ALL,REGION=caudate ../_h/step_1.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_h"
REGION="${REGION:-caudate}"
POPULATION="${POPULATION:-AA}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

python3 "${H}/01_define_scz_loci.py" --region "${REGION}" --population "${POPULATION}"
python3 "${H}/02_link_vmrs_to_loci.py" --region "${REGION}" --population "${POPULATION}"

log_message "**** Job ends ****"
