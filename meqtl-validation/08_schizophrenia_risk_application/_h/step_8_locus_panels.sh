#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=scz_locus_fig
#SBATCH --output=logs/scz_locus_fig.%j.log

# Prepare tidy locus tables + manuscript locus panels.
# Submit from meqtl-validation/08_schizophrenia_risk_application/_m/:
#   sbatch ../_h/step_8_locus_panels.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
python3 "${H}/21_prepare_locus_panel_data.py"

# ggplot2 / patchwork live in rnaseq
/projects/p32505/opt/envs/rnaseq/bin/Rscript "${H}/22_plot_locus_panels.R"

log_message "**** Job ends ****"
