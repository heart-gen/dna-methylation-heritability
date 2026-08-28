#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=array_universe
#SBATCH --output=logs/array_universe.%A.log
#
# inputs/supportfiles step 1: build an Illumina array probe universe in hg38.
#
# Produces the reference annotation that lets Module 01 quantify WGBS coverage
# outside array-accessible CpGs (AGENTS.md 2.2). Coordinates only -- no array
# intensities, samples, or methylation values are involved.
#
# 450K is the required comparator for Figure 1; EPIC is the stricter supplement.
# The EPIC annotation package must be installed in the epigenomics env:
#   BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#
# Usage, from the support-file _m directory:
#   cd inputs/supportfiles/_m && mkdir -p logs
#   PLATFORM=EPIC sbatch ../_h/step_1_build_array_universe.sh
#
# Prints the annotation_asset_manifest.tsv row to add; the manifest is edited
# deliberately, with the checksum verified, rather than appended by the job.

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${PLATFORM:?set PLATFORM=450K|EPIC}"

# liftOver and the hg19->hg38 chain are the only external dependencies.
require_exec "/projects/p32505/opt/envs/genomics/bin/liftOver"
require_file "$REPO_DIR/meqtl-validation/03_external_meqtl_validation/_m/support/hg19ToHg38.over.chain.gz"

log_job_info
log_message "**** Building ${PLATFORM} probe universe ****"

run_r "$REPO_DIR/inputs/supportfiles/_h/01_build_array_universe.R" \
    --platform "$PLATFORM"

log_message "**** Job ends ****"
