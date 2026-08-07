#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_phase4_sum
#SBATCH --output=logs/phase4_summary.%j.log

# Submit from meqtl-validation/04_cross_region_sharing/_m/ after step_1 and donor-group step:
#   sbatch --dependency=afterok:<cross>:<donor> ../_h/step_2.sh

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/04_cross_region_sharing/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Summarizing Phase 4 claims"
python3 "${H}/02_summarize_phase4.py"

log_message "**** Job ends ****"
