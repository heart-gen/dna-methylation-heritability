#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=dnam_scmd_validate
#SBATCH --output=logs/dnam_scmd_validate.%j.log

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/inputs/cell_proportions/_h"
mkdir -p logs

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
ENV_PREFIX="${CELL_DECONV_ENV:-${ROOT}/inputs/cell_proportions/_m/conda_env}"
if [[ ! -x "${ENV_PREFIX}/bin/Rscript" ]]; then
  echo "Missing ${ENV_PREFIX}; create it from conda_environments/cell_deconvolution.yml" >&2
  exit 1
fi
conda activate "${ENV_PREFIX}"
Rscript "${H}/07.validate_dnam_cell_prop.R"
Rscript "${H}/08.plot_dnam_cell_prop.R"
