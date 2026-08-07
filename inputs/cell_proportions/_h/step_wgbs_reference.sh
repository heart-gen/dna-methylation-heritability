#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --job-name=scmd_wgbs_ref
#SBATCH --output=logs/scmd_wgbs_ref.%j.log

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
ENV_PREFIX="${CELL_DECONV_ENV:-${ROOT}/inputs/cell_proportions/_m/conda_env}"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate "${ENV_PREFIX}"
Rscript "${ROOT}/inputs/cell_proportions/_h/04.prepare_scmd_markers.R" WGBS
