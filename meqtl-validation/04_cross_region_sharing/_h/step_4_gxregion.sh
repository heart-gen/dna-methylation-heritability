#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=phase4_gxregion
#SBATCH --output=logs/gxregion.%j.log

# Phase 4 Experiment 3: shared-donor genotype × region architecture screen.
# Submit from meqtl-validation/04_cross_region_sharing/_m/:
#   mkdir -p logs && sbatch ../_h/step_4_gxregion.sh
# Optional: MAX_PAIRS=1500 sbatch --export=ALL,MAX_PAIRS ../_h/step_4_gxregion.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/04_cross_region_sharing/_h"
MAX_PAIRS="${MAX_PAIRS:-3000}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Shared-donor G×region (max_pairs=${MAX_PAIRS})"
python3 "${H}/06_shared_donor_gxregion_meqtl.py" --max-pairs "${MAX_PAIRS}"

log_message "**** Job ends ****"
