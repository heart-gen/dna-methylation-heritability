#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=vmr_m6d
#SBATCH --output=logs/vmr_m6d.%A_%a.log

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h"
REGIONS=(caudate dlpfc hippocampus)
REGION="${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}"
CIS="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/covariate_sensitivity/tensorqtl/M6d/qc/lead_snp_per_cpg.tsv.gz"
OUT="${ROOT}/meqtl-validation/02_vmr_meqtl_burden/_m/M6d/${REGION}"
mkdir -p logs

if [[ ! -f "${CIS}" ]]; then
  echo "No gated M6d result for ${REGION}; skipping burden analysis."
  exit 0
fi
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
python3 "${H}/01_aggregate_vmr_burden.py" --region "${REGION}" --population AA \
  --cis-qtl "${CIS}" --outdir "${OUT}" --fdr 0.05
python3 "${H}/02_fit_burden_models.py" --region "${REGION}" \
  --burden-tsv "${OUT}/vmr_meqtl_burden.tsv.gz"
