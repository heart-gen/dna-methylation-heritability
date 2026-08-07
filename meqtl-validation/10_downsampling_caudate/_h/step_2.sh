#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-29%8
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=ds_caud_tqtl
#SBATCH --output=logs/tensorqtl_ds.%A_%a.log

# GPU array: official TensorQTL cis permutation FDR per downsample replicate.
# Requires step_1. Submit from meqtl-validation/10_downsampling_caudate/_m/:
#   sbatch --dependency=afterok:<prep> ../_h/step_2.sh
# Smoke: sbatch --array=0 --export=ALL,MAX_REPS=1 ../_h/step_2.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}; array: ${SLURM_ARRAY_TASK_ID:-N/A}; CUDA=${CUDA_VISIBLE_DEVICES:-unset}"
module purge
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H10="${ROOT}/meqtl-validation/10_downsampling_caudate/_h"
M10="${ROOT}/meqtl-validation/10_downsampling_caudate/_m"
H01="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
GENO="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/caudate/_m/genotypes/meqtl_AA"
DEVICE="${DEVICE:-gpu}"
CHUNK_SIZE="${CHUNK_SIZE:-chr}"
MAX_REPS="${MAX_REPS:-0}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
python -c "import tensorqtl" >/dev/null 2>&1 \
  || { log_message "ERROR: tensorqtl missing"; exit 1; }

MANIFEST="${M10}/downsample_replicate_manifest.tsv"
if [[ ! -f "${MANIFEST}" ]]; then
  log_message "ERROR: missing ${MANIFEST}; run step_1 first"
  exit 1
fi

# Array index 0 -> replicate 1
REP=$(( SLURM_ARRAY_TASK_ID + 1 ))
if [[ "${MAX_REPS}" != "0" && "${REP}" -gt "${MAX_REPS}" ]]; then
  log_message "REP ${REP} > MAX_REPS ${MAX_REPS}; exiting"
  exit 0
fi

# Pull paths from manifest (awk: replicate is col1)
PHENO=$(awk -F'\t' -v r="${REP}" 'NR==1{for(i=1;i<=NF;i++) if($i=="phenotype_bed") p=i; next} $1==r{print $p; exit}' "${MANIFEST}")
COV=$(awk -F'\t' -v r="${REP}" 'NR==1{for(i=1;i<=NF;i++) if($i=="covariates") p=i; next} $1==r{print $p; exit}' "${MANIFEST}")
OUTDIR=$(awk -F'\t' -v r="${REP}" 'NR==1{for(i=1;i<=NF;i++) if($i=="outdir") p=i; next} $1==r{print $p; exit}' "${MANIFEST}")
PREFIX=$(awk -F'\t' -v r="${REP}" 'NR==1{for(i=1;i<=NF;i++) if($i=="prefix") p=i; next} $1==r{print $p; exit}' "${MANIFEST}")

if [[ -z "${PHENO}" || -z "${COV}" || -z "${OUTDIR}" ]]; then
  log_message "ERROR: could not resolve manifest row for replicate ${REP}"
  exit 1
fi
if [[ ! -f "${PHENO}" ]]; then
  log_message "ERROR: missing phenotype BED ${PHENO}"
  exit 1
fi
if [[ ! -f "${COV}" ]]; then
  log_message "ERROR: missing covariates ${COV}"
  exit 1
fi
if [[ ! -f "${GENO}.pgen" ]]; then
  log_message "ERROR: missing genotypes ${GENO}.pgen"
  exit 1
fi

mkdir -p "${OUTDIR}"
log_message "TensorQTL cis map rep=${REP} prefix=${PREFIX} device=${DEVICE}"
python "${H01}/04_tensorqtl_map.py" \
  --region caudate \
  --mode cis \
  --phenotype-bed "${PHENO}" \
  --covariates "${COV}" \
  --genotype-prefix "${GENO}" \
  --outdir "${OUTDIR}" \
  --prefix "${PREFIX}" \
  --window 500000 \
  --maf 0.05 \
  --seed 20260805 \
  --device "${DEVICE}" \
  --chunk-size "${CHUNK_SIZE}"

CIS="${OUTDIR}/${PREFIX}.cis_qtl.txt.gz"
if [[ ! -f "${CIS}" ]]; then
  log_message "ERROR: missing cis output ${CIS}"
  exit 1
fi

log_message "QC summarize rep=${REP}"
python "${H01}/05_qc_summarize.py" \
  --region caudate \
  --cis-qtl "${CIS}" \
  --outdir "${OUTDIR}/qc"

log_message "**** Job ends ****"
