#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=dnam_scmd_ref
#SBATCH --output=logs/dnam_scmd_ref.%j.log

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

bash "${H}/03.download_scmd_reference.sh"
Rscript "${H}/03a.install_annotation.R"
if [[ "${INSTALL_SCMD:-1}" == "1" ]]; then
  if Rscript "${H}/03b.install_scmd.R"; then
    printf 'status\tnote\ninstalled\texact upstream scMD available\n' \
      > "${ROOT}/inputs/cell_proportions/_m/reference/scmd_install_status.tsv"
  else
    printf 'status\tnote\nfailed\tsee job log; recorded fallback remains enabled\n' \
      > "${ROOT}/inputs/cell_proportions/_m/reference/scmd_install_status.tsv"
    if [[ "${REQUIRE_UPSTREAM_SCMD:-0}" == "1" ]]; then exit 1; fi
  fi
fi
Rscript "${H}/04.prepare_scmd_markers.R"
