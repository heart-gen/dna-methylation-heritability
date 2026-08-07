#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=scz_caud_down
#SBATCH --output=logs/scz_caud_downsample.%j.log

# Tier A: caudate donor-downsample targeted SCZ meQTL.
# Submit from meqtl-validation/08_schizophrenia_risk_application/_m/:
#   sbatch ../_h/step_4.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_h"
N_REPS="${N_REPS:-30}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Build downsample sample lists (n_reps=${N_REPS})"
python3 "${H}/08_make_caudate_downsample_lists.py" --n-reps "${N_REPS}"

log_message "Run downsampled targeted meQTL"
python3 "${H}/09_run_caudate_downsample_meqtl.py"

log_message "Summarize vs DLPFC/hippocampus"
python3 "${H}/10_summarize_caudate_downsample.py"

log_message "**** Job ends ****"
