#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=celltype_line
#SBATCH --output=logs/celltype_line.%j.log

# Cell-type LINE/L1 / repressive-compartment sensitivity.
# Submit from meqtl-validation/11_celltype_compartment_sensitivity/_m/:
#   sbatch ../_h/step_1.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/11_celltype_compartment_sensitivity/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

python3 "${H}/01_audit_celltype_overlap.py"
python3 "${H}/02_build_vmr_cell_metrics.py"
python3 "${H}/03_run_celltype_sensitivity.py"

log_message "**** Job ends ****"
