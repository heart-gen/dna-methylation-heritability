#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_external_harm
#SBATCH --output=logs/external_harmonize.%j.log

# Submit from meqtl-validation/03_external_meqtl_validation/_m:
#   sbatch ../_h/step_harmonize.sh

set -euo pipefail

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/03_external_meqtl_validation/_h"
cd "${ROOT}/meqtl-validation/03_external_meqtl_validation/_m"
mkdir -p logs

echo "$(date '+%Y-%m-%d %H:%M:%S') - **** Job starts ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
export PYTHONPATH="/home/owb0346/.local/lib/python3.11/site-packages:${PYTHONPATH:-}"

python3 "${H}/02_download_and_record.py"
python3 "${H}/03_harmonize_external_meqtls.py" \
  --resources jaffe_dlpfc_450k_meqtl schulz_hippocampus_array_meqtl brainseq_wgbs_meqtl_scz_subset \
  --fdr 0.05

echo "$(date '+%Y-%m-%d %H:%M:%S') - **** Job ends ****"
