#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_cpg_tensorqtl
#SBATCH --output=logs/cpg_tensorqtl.%A_%a.log

# GPU: TensorQTL cis-meQTL mapping
# Requires step_2 phenotypes and step_3 genotypes. Submit from _m/:
#   sbatch ../_h/step_4.sh
# Fallback CPU (if GPU unavailable):
#   sbatch --partition=genomics --account=b1042 --gres= --time=2-00:00:00 \
#     --export=ALL,DEVICE=cpu ../_h/step_4.sh

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Node name: ${SLURM_NODENAME:-N/A}"
echo "Array task: ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"

module purge
module list
mkdir -p logs

REGIONS=(caudate dlpfc hippocampus)
if [[ -z "${REGION:-}" ]]; then
  REGION="${REGIONS[${SLURM_ARRAY_TASK_ID}]}"
fi
POPULATION="${POPULATION:-AA}"
DEVICE="${DEVICE:-gpu}"
CHUNK_SIZE="${CHUNK_SIZE:-chr}"

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
# Submit from _m/; do not resolve $0 (SLURM copies the batch script to spool).
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
M="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m"
if [[ "${POPULATION}" == "AA" ]]; then
  PREP="${M}/prepared"
  OUTDIR="${M}/tensorqtl"
  PREFIX="cpg_meqtl_${REGION}"
else
  PREP="${M}/prepared/${POPULATION}"
  OUTDIR="${M}/tensorqtl/${POPULATION}"
  PREFIX="cpg_meqtl_${REGION}_${POPULATION}"
fi
PHENO_BED="${PREP}/cpg_phenotypes.all_autosomes.bed"
PHENO_GZ="${PREP}/cpg_phenotypes.all_autosomes.bed.gz"
GENO_PREFIX="${M}/genotypes/meqtl_${POPULATION}"

log_message "**** Loading conda environment ****"
# shellcheck disable=SC1091
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
python -c "import tensorqtl" >/dev/null 2>&1 \
  || { log_message "ERROR: tensorqtl missing in genomics env"; exit 1; }

# Ensure bgzipped phenotype BED (tensorqtl accepts .bed, but gz is preferred)
if [[ ! -f "${PHENO_GZ}" ]]; then
  if [[ ! -f "${PHENO_BED}" ]]; then
    log_message "ERROR: missing phenotype BED: ${PHENO_BED}"
    exit 1
  fi
  log_message "bgzip + tabix phenotype BED for ${REGION} (${POPULATION})"
  bgzip -c "${PHENO_BED}" > "${PHENO_GZ}"
  tabix -f -p bed "${PHENO_GZ}"
fi
PHENO="${PHENO_GZ}"

if [[ ! -f "${GENO_PREFIX}.pgen" ]]; then
  log_message "ERROR: missing filtered genotypes: ${GENO_PREFIX}.pgen"
  exit 1
fi
if [[ ! -f "${PREP}/covariates.txt" ]]; then
  log_message "ERROR: missing covariates: ${PREP}/covariates.txt"
  exit 1
fi

mkdir -p "${OUTDIR}"
log_message "Running TensorQTL cis-meQTL for ${REGION} (${POPULATION}; device=${DEVICE}, chunk=${CHUNK_SIZE})"
python "${H}/04_tensorqtl_map.py" \
  --region "${REGION}" \
  --mode cis \
  --phenotype-bed "${PHENO}" \
  --covariates "${PREP}/covariates.txt" \
  --genotype-prefix "${GENO_PREFIX}" \
  --outdir "${OUTDIR}" \
  --prefix "${PREFIX}" \
  --window 500000 \
  --maf 0.05 \
  --device "${DEVICE}" \
  --chunk-size "${CHUNK_SIZE}"

conda deactivate
log_message "**** Job ends ****"
