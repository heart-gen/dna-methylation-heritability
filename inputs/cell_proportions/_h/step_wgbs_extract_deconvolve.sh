#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-1
#SBATCH --job-name=scmd_wgbs_sens
#SBATCH --output=logs/scmd_wgbs_sens.%A_%a.log

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/inputs/cell_proportions/_h"
REGIONS=(dlpfc hippocampus)
REGION="${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}"
ENV_PREFIX="${CELL_DECONV_ENV:-${ROOT}/inputs/cell_proportions/_m/conda_env}"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate "${ENV_PREFIX}"
Rscript "${H}/05.extract_wgbs_markers.R" "${REGION}" WGBS
Rscript "${H}/06.run_scmd.R" "${REGION}" WGBS
