#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=libd_pc5_cov
#SBATCH --output=logs/libd_pc5_cov.%j.log

# Sensitivity: reuse RPKM prepared expression; rebuild covariates with --max-pcs 5
# Does not overwrite genes/ or genes_rpkm/standard.

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_h"
REGION="${REGION:-caudate}"
PREP="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_m/${REGION}/genes_rpkm/prepared"
OUT="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_m/${REGION}/genes_rpkm_pc5"
STD="${OUT}/standard"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/rnaseq

[[ -f "${PREP}/normalized_expression.tsv.gz" ]] || { log_message "ERROR: missing RPKM prepared expression"; exit 1; }
mkdir -p "${STD}"
# Reuse phenotype BED from RPKM prep via symlink
mkdir -p "${OUT}"
ln -sfn "${PREP}" "${OUT}/prepared"

Rscript "${H}/02_make_covariates.R" \
  --prepared-dir "${PREP}" \
  --outdir "${STD}" \
  --max-pcs 5

log_message "**** Job ends ****"
