#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=dnam_scmd_extract
#SBATCH --output=logs/dnam_scmd_extract.%A_%a.log

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/inputs/cell_proportions/_h"
REGIONS=(caudate dlpfc hippocampus)
REGION="${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}"
mkdir -p logs

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
ENV_PREFIX="${CELL_DECONV_ENV:-${ROOT}/inputs/cell_proportions/_m/conda_env}"
if [[ ! -x "${ENV_PREFIX}/bin/Rscript" ]]; then
  echo "Missing ${ENV_PREFIX}; create it from conda_environments/cell_deconvolution.yml" >&2
  exit 1
fi
conda activate "${ENV_PREFIX}"
Rscript "${H}/05.extract_wgbs_markers.R" "${REGION}"
