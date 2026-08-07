#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=scz_l3_eqtl
#SBATCH --output=logs/scz_l3_eqtl.%j.log

# Compatibility wrapper: Level-3 eQTL tests now live in step_6b_level3_eqtl_tests.sh
# Genome-wide LIBD eQTL mapping: meqtl-validation/09_libd_eqtl_mapping/

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
exec bash "${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_h/step_6b_level3_eqtl_tests.sh"
