#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=phase4_caud_down
#SBATCH --output=logs/caudate_downsample.%j.log

# Phase 4 Experiment 3: caudate N-matched lead-SNP retention downsample.
# Submit from meqtl-validation/04_cross_region_sharing/_m/:
#   mkdir -p logs && sbatch ../_h/step_3_downsample.sh
# Optional: N_REPS=10 MAX_REPS=10 sbatch --export=ALL,N_REPS,MAX_REPS ../_h/step_3_downsample.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/04_cross_region_sharing/_h"
N_REPS="${N_REPS:-30}"
MAX_REPS="${MAX_REPS:-0}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Build downsample sample lists (n_reps=${N_REPS})"
python3 "${H}/03_make_caudate_downsample_lists.py" --n-reps "${N_REPS}"

log_message "Run N-matched lead-SNP retention meQTL"
if [[ "${MAX_REPS}" != "0" ]]; then
  python3 "${H}/04_run_caudate_downsample_meqtl.py" --max-reps "${MAX_REPS}"
else
  python3 "${H}/04_run_caudate_downsample_meqtl.py"
fi

log_message "Summarize vs DLPFC/hippocampus"
python3 "${H}/05_summarize_caudate_downsample.py"

log_message "**** Job ends ****"
