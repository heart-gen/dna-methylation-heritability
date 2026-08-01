#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_covsens
#SBATCH --output=logs/covsens_tensorqtl.%A_%a.log

# TensorQTL for covariate-sensitivity models (does NOT overwrite primary tensorqtl/).
# Default REGION=caudate pilot models: M1 M2 M3a M3b M3c
# Submit from _m/ after step_8:
#   sbatch --export=ALL,REGION=caudate ../_h/step_9_covar_models.sh
# CPU fallback:
#   sbatch --partition=genomics --account=b1042 --gres= --time=2-00:00:00 \
#     --export=ALL,REGION=caudate,DEVICE=cpu ../_h/step_9_covar_models.sh
# Custom model list:
#   sbatch --export=ALL,REGION=dlpfc,MODELS="M1 M3b" --array=0-1 ../_h/step_9_covar_models.sh

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

REGION="${REGION:-caudate}"
DEVICE="${DEVICE:-gpu}"
CHUNK_SIZE="${CHUNK_SIZE:-chr}"
# shellcheck disable=SC2206
MODELS=(${MODELS:-M1 M2 M3a M3b M3c})
MODEL="${MODELS[${SLURM_ARRAY_TASK_ID}]}"

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
M="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m"
PREP="${M}/prepared"
SENS="${M}/covariate_sensitivity"
COV="${SENS}/covariates_${MODEL}.txt"
OUT="${SENS}/tensorqtl/${MODEL}"
PHENO_GZ="${PREP}/cpg_phenotypes.all_autosomes.bed.gz"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

if [[ ! -f "${COV}" ]]; then
  log_message "ERROR: missing ${COV}; run step_8_latent.sh first"
  exit 1
fi
if [[ ! -f "${PHENO_GZ}" ]]; then
  log_message "ERROR: missing phenotype BED ${PHENO_GZ}"
  exit 1
fi
if [[ ! -f "${M}/genotypes/meqtl_AA.pgen" ]]; then
  log_message "ERROR: missing genotypes"
  exit 1
fi

mkdir -p "${OUT}/qc"

log_message "TensorQTL cis ${REGION} model=${MODEL} device=${DEVICE}"
python3 "${H}/04_tensorqtl_map.py" \
  --region "${REGION}" \
  --mode cis \
  --phenotype-bed "${PHENO_GZ}" \
  --covariates "${COV}" \
  --genotype-prefix "${M}/genotypes/meqtl_AA" \
  --outdir "${OUT}" \
  --prefix "cpg_meqtl_${REGION}_${MODEL}" \
  --window 500000 \
  --maf 0.05 \
  --device "${DEVICE}" \
  --chunk-size "${CHUNK_SIZE}"

CIS="${OUT}/cpg_meqtl_${REGION}_${MODEL}.cis_qtl.txt.gz"
log_message "QC summarize ${MODEL}"
python3 "${H}/05_qc_summarize.py" \
  --region "${REGION}" \
  --cis-qtl "${CIS}" \
  --outdir "${OUT}/qc"

log_message "Calibration plots ${MODEL}"
python3 "${H}/08_calibration_plots.py" \
  --region "${REGION}" \
  --cis-qtl "${CIS}" \
  --outdir "${OUT}/qc/calibration"

log_message "**** Job ends ****"
