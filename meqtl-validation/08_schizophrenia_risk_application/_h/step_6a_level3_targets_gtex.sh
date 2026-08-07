#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=scz_l3_gtex
#SBATCH --output=logs/scz_l3_gtex.%j.log

# Level 3: gene targets + GTEx lookup
# Submit from meqtl-validation/08_schizophrenia_risk_application/_m/:
#   sbatch ../_h/step_6a_level3_targets_gtex.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

python3 "${H}/12_build_level3_gene_targets.py"
python3 "${H}/13_gtex_level3_lookup.py"

log_message "**** Job ends ****"
