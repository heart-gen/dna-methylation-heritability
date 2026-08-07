#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_repeat_robust
#SBATCH --output=logs/repeat_robust.%j.log

# Submit from meqtl-validation/07_repeat_mappability_sensitivity/_m/:
#   mkdir -p logs && sbatch ../_h/step_1.sh

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Running consolidated repeat/mappability robustness analyses"
python3 "${H}/03_run_robustness_analyses.py"

log_message "**** Job ends ****"
