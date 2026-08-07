#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=ds_caud_sum
#SBATCH --output=logs/summarize_downsample.%j.log

# CPU: aggregate official TensorQTL downsample vs DLPFC/hippocampus.
# Submit from meqtl-validation/10_downsampling_caudate/_m/:
#   sbatch --dependency=afterok:<array> ../_h/step_3.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/10_downsampling_caudate/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Summarizing TensorQTL downsample replicates"
python3 "${H}/04_summarize_tensorqtl_downsample.py"

log_message "**** Job ends ****"
