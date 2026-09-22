#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=scz-negative-controls
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/scz_negative_controls.%A.log
#
# 09 stages 17 and 18 -- the trait-general qualification of the Module 09 axis
# result (AGENTS.md 7.8), over the 114-trait harmonized GWAS collection.
#
#   stage 17  locus -> VMR -> axis contrast for every trait, read by category
#   stage 18  which Module 04 annotation separates the depleted traits
#
# Both are MODULE-LEVEL: they read the accepted per-region runs and write
# _m/combined/. They mint no run ID and change no decision.
#
# Why this script exists. Both stages were hand-run before, and stage 17 holds
# the lead sets for all 117 extracted tags in memory while it fits; on a login
# node it is OOM-killed part way through the trait loop. It needs a batch
# allocation, and a launcher is also what makes the re-run reproducible.
#
# Prerequisite: 17a_extract_gwas_leads.sh must have completed, which stage 17
# checks by the sig/.complete marker.
#
# ---------------------------------------------------------------------------
# ACCEPTED vs UNACCEPTED output
#
# Neither stage takes a scientific switch. `--allow-unlocked` relaxes TWO
# things at once -- assert_locked() on config/gwas_negative_controls.yml, and
# require_accepted_upstream() on the Module 09 runs -- and it is the ONLY thing
# that appends `-UNACCEPTED` to the output file names.
#
# config/gwas_negative_controls.yml was PI-locked on 2026-09-20 and the three
# Module 09 runs scz-AA-{caudate,dlpfc,hippocampus}-20260918 were accepted on
# 2026-09-19, so both assertions now pass on their own. This script therefore
# NEVER passes the flag: the outputs are citable or the stage stops. Set
# ALLOW_UNLOCKED=1 only to reproduce the superseded uncitable tables.
#
# Usage, from the module's _m directory:
#   cd 09_schizophrenia_risk_application/_m && mkdir -p logs
#   sbatch ../_h/step_9_negative_controls.sh
#   COHORT=AA sbatch ../_h/step_9_negative_controls.sh

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

COHORT="${COHORT:-AA}"
HERE="$REPO_DIR/09_schizophrenia_risk_application/_h"

EXTRA=()
if [ "${ALLOW_UNLOCKED:-0}" = "1" ]; then
    EXTRA+=(--allow-unlocked)
    log_message "WARNING: --allow-unlocked set; output will carry -UNACCEPTED and may not be cited"
fi

log_job_info
log_message "**** Module 09 negative-control stages, cohort=${COHORT} ****"

log_message "[17] negative-control traits"
run_r "$HERE/17_negative_control_traits.R" --cohort "$COHORT" "${EXTRA[@]}"

log_message "[18] locus architecture"
run_r "$HERE/18_locus_architecture.R" --cohort "$COHORT" "${EXTRA[@]}"

log_message "**** Job ends ****"
