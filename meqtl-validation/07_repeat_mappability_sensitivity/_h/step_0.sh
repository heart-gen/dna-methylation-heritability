#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_ann_assets
#SBATCH --output=logs/annotation_assets.%j.log

# Submit from meqtl-validation/07_repeat_mappability_sensitivity/_m:
#   mkdir -p logs && sbatch ../_h/step_0_assets.sh

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_h"
cd "${ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_m"
mkdir -p logs

echo "$(date '+%Y-%m-%d %H:%M:%S') - **** Job starts ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

# Assets already staged by prior run; --force not set so downloads are skipped.
python3 "${H}/00_stage_annotation_assets.py"
python3 "${H}/02_annotate_vmr_technical.py"
python3 "${H}/01_build_robustness_table.py"

echo "$(date '+%Y-%m-%d %H:%M:%S') - **** Job ends ****"
