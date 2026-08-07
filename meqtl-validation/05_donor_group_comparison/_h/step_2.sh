#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_donor_depth
#SBATCH --output=logs/donor_group_depth.%j.log

# Phase 4b Experiment 2 depth: AA–EA effect concordance + MAF/LD-matched discovery.
# Submit from meqtl-validation/05_donor_group_comparison/_m/:
#   mkdir -p logs && sbatch ../_h/step_2.sh

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/05_donor_group_comparison/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "AA–EA CpG lead-effect concordance"
python3 "${H}/02_aa_ea_effect_concordance.py" --region all

log_message "MAF / cis-SNP-density matched discovery contrasts"
python3 "${H}/03_maf_ld_matched_discovery.py" --region all

log_message "**** Job ends ****"
