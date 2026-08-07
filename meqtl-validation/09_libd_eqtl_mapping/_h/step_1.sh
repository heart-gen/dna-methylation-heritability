#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=libd_eqtl_prep
#SBATCH --output=logs/libd_eqtl_prep.%j.log

# AA-only LIBD expression prep + covariates + phenotype BED
# Submit from meqtl-validation/09_libd_eqtl_mapping/_m/:
#   sbatch --export=ALL,REGION=caudate ../_h/step_1.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_h"
REGION="${REGION:-caudate}"
OUT="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_m/${REGION}/genes"
PREP="${OUT}/prepared"
STD="${OUT}/standard"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/rnaseq

mkdir -p "${PREP}" "${STD}"
Rscript "${H}/01_prepare_expression.R" \
  --region "${REGION}" \
  --outdir "${PREP}" \
  --min-age 13

Rscript "${H}/02_make_covariates.R" \
  --prepared-dir "${PREP}" \
  --outdir "${STD}"

conda activate /projects/p32505/opt/envs/genomics
python3 "${H}/03_make_tensorqtl_bed.py" \
  --prepared-dir "${PREP}" \
  --outdir "${PREP}" \
  --prefix genes

log_message "**** Job ends ****"
