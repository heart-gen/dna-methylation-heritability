#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_external_init
#SBATCH --output=logs/external_init.%j.log

# Submit from meqtl-validation/03_external_meqtl_validation/_m/:
#   mkdir -p logs && sbatch ../_h/step_1.sh

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/03_external_meqtl_validation/_h"

module purge
module list

log_message "Initializing external meQTL validation workspace"
python3 "${H}/00_init_workspace.py"

log_message "**** Job ends ****"
