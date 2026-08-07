#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_m6d
#SBATCH --output=logs/m6d_tensorqtl.%A_%a.log

# Gated DNAm-composition sensitivity. This never writes into primary M3a paths.
set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
REGIONS=(caudate dlpfc hippocampus)
REGION="${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}"
DEVICE="${DEVICE:-gpu}"
M="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m"
PREP="${M}/prepared"
SENS="${M}/covariate_sensitivity"
COV="${SENS}/covariates_M6d.txt"
SAMPLE_LIST="${SENS}/dnam_cell_pca/m6d_sample_ids.txt"
mkdir -p logs

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
python3 "${H}/10_prepare_covariate_models.py" --region "${REGION}"
if [[ ! -f "${COV}" || ! -f "${SAMPLE_LIST}" ]]; then
  echo "M6d was not created for ${REGION}; the integration gate failed or required inputs are absent."
  exit 0
fi

for MODEL in M3a_dnam_matched M6d; do
  MODEL_COV="${SENS}/covariates_${MODEL}.txt"
  OUT="${SENS}/tensorqtl/${MODEL}"
  mkdir -p "${OUT}/qc"
  python3 "${H}/04_tensorqtl_map.py" \
    --region "${REGION}" --mode cis \
    --phenotype-bed "${PREP}/cpg_phenotypes.all_autosomes.bed.gz" \
    --covariates "${MODEL_COV}" --genotype-prefix "${M}/genotypes/meqtl_AA" \
    --sample-list "${SAMPLE_LIST}" \
    --outdir "${OUT}" --prefix "cpg_meqtl_${REGION}_${MODEL}" \
    --window 500000 --maf 0.05 --device "${DEVICE}" --chunk-size chr
  CIS="${OUT}/cpg_meqtl_${REGION}_${MODEL}.cis_qtl.txt.gz"
  python3 "${H}/05_qc_summarize.py" --region "${REGION}" --cis-qtl "${CIS}" --outdir "${OUT}/qc"
  python3 "${H}/08_calibration_plots.py" --region "${REGION}" --cis-qtl "${CIS}" --outdir "${OUT}/qc/calibration"
done
python3 "${H}/11_compare_covariate_models.py" --region "${REGION}" \
  --models M3a_dnam_matched M6d --baseline-model M3a_dnam_matched \
  --output-name m6d_vs_m3a_matched_summary.tsv
