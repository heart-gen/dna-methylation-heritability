#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_external_enrich
#SBATCH --output=logs/external_enrich.%j.log

# Submit from meqtl-validation/03_external_meqtl_validation/_m/:
#   mkdir -p logs && sbatch ../_h/step_2.sh
# Or run locally with genomics conda after Phase 2 burden tables exist.

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/03_external_meqtl_validation/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Running Phase 3 external enrichment for all resources × regions"
python3 "${H}/04_run_and_summarize.py"

log_message "**** Job ends ****"
