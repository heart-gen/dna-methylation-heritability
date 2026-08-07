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
#SBATCH --job-name=ea_m3a_tensorqtl
#SBATCH --output=logs/ea_m3a_tensorqtl.%j.log

# EA-M3a TensorQTL for one region (default caudate).
# Submit from _m/ after step_ea_m3a_caudate.sh:
#   sbatch --dependency=afterok:<latent> ../_h/step_ea_m3a_tensorqtl.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
REGION="${REGION:-caudate}"
POPULATION=EA
DEVICE="${DEVICE:-gpu}"
CHUNK_SIZE="${CHUNK_SIZE:-chr}"

M="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m"
PREP="${M}/prepared/${POPULATION}"
PHENO="${PREP}/cpg_phenotypes.all_autosomes.bed.gz"
COV="${PREP}/covariates_M3a.txt"
GENO="${M}/genotypes/meqtl_${POPULATION}"
OUTDIR="${M}/tensorqtl/${POPULATION}/M3a"
PREFIX="cpg_meqtl_${REGION}_${POPULATION}_M3a"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

[[ -f "${PHENO}" ]] || { log_message "ERROR: missing ${PHENO}"; exit 1; }
[[ -f "${COV}" ]] || { log_message "ERROR: missing ${COV}"; exit 1; }
[[ -f "${GENO}.pgen" ]] || { log_message "ERROR: missing ${GENO}.pgen"; exit 1; }

mkdir -p "${OUTDIR}"
log_message "EA-M3a TensorQTL ${REGION} (device=${DEVICE})"
head -1 "${COV}"

python "${H}/04_tensorqtl_map.py" \
  --region "${REGION}" \
  --mode cis \
  --phenotype-bed "${PHENO}" \
  --covariates "${COV}" \
  --genotype-prefix "${GENO}" \
  --outdir "${OUTDIR}" \
  --prefix "${PREFIX}" \
  --window 500000 \
  --maf 0.05 \
  --device "${DEVICE}" \
  --chunk-size "${CHUNK_SIZE}"

python "${H}/05_qc_summarize.py" \
  --region "${REGION}" \
  --cis-qtl "${OUTDIR}/${PREFIX}.cis_qtl.txt.gz"

log_message "**** Job ends ****"
