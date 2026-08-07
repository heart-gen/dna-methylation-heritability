#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=ds_caud_prep
#SBATCH --output=logs/prepare_downsample.%j.log

# CPU: prepare Exp3 sample lists + per-rep covariates + phenotype BEDs.
# Submit from meqtl-validation/10_downsampling_caudate/_m/:
#   mkdir -p logs && sbatch ../_h/step_1.sh
# Optional: MAX_REPS=1 sbatch --export=ALL,MAX_REPS ../_h/step_1.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/10_downsampling_caudate/_h"
MAX_REPS="${MAX_REPS:-0}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Preparing downsample inputs (MAX_REPS=${MAX_REPS})"
if [[ "${MAX_REPS}" != "0" ]]; then
  python3 "${H}/01_prepare_downsample_inputs.py" --max-reps "${MAX_REPS}"
else
  python3 "${H}/01_prepare_downsample_inputs.py"
fi

log_message "**** Job ends ****"
