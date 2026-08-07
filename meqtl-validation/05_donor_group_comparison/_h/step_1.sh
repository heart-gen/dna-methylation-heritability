#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_donor_group
#SBATCH --output=logs/donor_group.%j.log

# Submit from meqtl-validation/05_donor_group_comparison/_m/:
#   mkdir -p logs && sbatch ../_h/step_1.sh
#
# EA stratified meQTL + burden are complete (M0; see ea_meqtl_readiness.tsv).
# This job runs predictability portability + readiness + coefficient comparison.
# Experiment 2 depth (concordance + MAF/LD matching): sbatch ../_h/step_2.sh

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/05_donor_group_comparison/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Donor-group predictability portability and EA readiness"
python3 "${H}/01_compare_donor_groups.py" --region all

log_message "**** Job ends ****"
