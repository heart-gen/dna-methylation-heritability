#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=libd_rpkm_tqtl
#SBATCH --output=logs/libd_rpkm_tqtl.%j.log

# Sensitivity TensorQTL for genes_rpkm/ (does not overwrite CPM primary)
# Submit from _m/ after step_1_rpkm_sensitivity.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
MAPPER="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h/04_tensorqtl_map.py"
REGION="${REGION:-caudate}"
OUT="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_m/${REGION}/genes_rpkm"
PREP="${OUT}/prepared"
STD="${OUT}/standard"
TQTL="${OUT}/tensorqtl"
PHENO="${PREP}/genes.expression.bed.gz"
COV="${STD}/covariates.txt"
GENO="${ROOT}/inputs/genotypes/TOPMed_LIBD.AA"
PREFIX="libd_aa_${REGION}_genes_rpkm"
DEVICE="${DEVICE:-gpu}"
CHUNK_SIZE="${CHUNK_SIZE:-chr}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

[[ -f "${PHENO}" ]] || { log_message "ERROR: missing ${PHENO}"; exit 1; }
[[ -f "${COV}" ]] || { log_message "ERROR: missing ${COV}"; exit 1; }

mkdir -p "${TQTL}"
log_message "RPKM sensitivity tensorQTL cis ${REGION} (device=${DEVICE})"
head -1 "${COV}"
python - <<PY
import pandas as pd
cov=pd.read_csv("${COV}",sep="\t",nrows=0)
pcs=[c for c in cov.columns if c.startswith("PC") and c[2:].isdigit()]
print("n_expr_pcs", len(pcs))
PY

python "${MAPPER}" \
  --region "${REGION}" \
  --mode cis \
  --phenotype-bed "${PHENO}" \
  --covariates "${COV}" \
  --genotype-prefix "${GENO}" \
  --outdir "${TQTL}" \
  --prefix "${PREFIX}" \
  --window 500000 \
  --maf 0.01 \
  --device "${DEVICE}" \
  --chunk-size "${CHUNK_SIZE}" \
  --seed 13131313 \
  --fdr 0.05

log_message "**** Job ends ****"
