#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=manuscript_figures
#SBATCH --output=logs/manuscript_figures.%A.log
#
# 11_integrated_manuscript_outputs step 1: assemble Figures 1 and 2.
#
# Consumes only accepted immutable upstream runs (AGENTS.md 7.9). Figure 1
# needs the Module 01 catalog runs AND the QC refresh runs that carry array
# coverage and genomic context; Figure 2 needs the Module 02 score runs.
#
# The builders resolve those upstream run IDs internally, so the only argument
# is the output run ID.
#
# Usage, from the module's _m directory:
#   cd 11_integrated_manuscript_outputs/_m && mkdir -p logs
#   RUN_ID=fig-all-20260826-a sbatch ../_h/step_1_figures.sh

# SLURM copies this script into a spool directory, so BASH_SOURCE does not
# point at the repository. Resolve the root from the submission directory.
_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${RUN_ID:?set RUN_ID}"

# Execute the run's own snapshot of _h/, not the live working tree. The submit
# driver copies _h/ and config/ into {run_dir}/code and exports V2_RUN_CODE, so
# an edit to _h/ after submission cannot change what this run built. The
# fallback keeps a hand-run build working, and such a run simply has no
# snapshot to cite.
HERE="${V2_RUN_CODE:-$REPO_DIR/11_integrated_manuscript_outputs/_h}"
mkdir -p "$REPO_DIR/11_integrated_manuscript_outputs/_m/runs/$RUN_ID"

log_job_info
log_message "**** Building manuscript figures into ${RUN_ID} ****"

# ---------------------------------------------------------------- both arms
# AA is the primary arm; all_individuals renders from the same builders as the
# sensitivity supplement. Figure 1 gets a second pass with EPIC as the stricter
# array comparator.
for COHORT in AA all_individuals; do
    log_message "Figure 1 (450K) + catalog-turnover supplement -- ${COHORT}"
    run_r "$HERE/01_figure1_catalog.R" --cohort "$COHORT" --run-id "$RUN_ID"

    log_message "Figure 1 (EPIC supplement) -- ${COHORT}"
    run_r "$HERE/01_figure1_catalog.R" --cohort "$COHORT" --run-id "$RUN_ID" \
        --platform EPIC

    log_message "Figure 2 + audit and denominator supplements -- ${COHORT}"
    run_r "$HERE/02_figure2_local_control.R" --cohort "$COHORT" --run-id "$RUN_ID"
done

# --------------------------------------------------------------- AA only
# Modules 04-10 have accepted runs for the AA arm only, so Figures 3-5 and the
# supplements below are NOT in the two-arm loop above. Running them for
# all_individuals would fail the acceptance gate, which is the correct
# behaviour but a confusing way to discover it.
log_message "Figure 3: repeat and repressive chromatin (Module 04) -- AA"
run_r "$HERE/07_figure3_repeat_repressive.R" --cohort AA --run-id "$RUN_ID"

log_message "Figure 4: meQTL burden and coupling (Modules 05, 07) -- AA"
run_r "$HERE/08_figure4_meqtl_coupling.R" --cohort AA --run-id "$RUN_ID"

log_message "Figure 5: trait-general GWAS architecture (Module 09 stages 17/18) -- AA"
run_r "$HERE/09_figure5_gwas_architecture.R" --cohort AA --run-id "$RUN_ID"

# Module 08 is AA-only and spans all three regions; the builder resolves the
# accepted crossregion run itself.
log_message "Region and donor-group generalization (Module 08) -- AA"
run_r "$HERE/06_figure_region_donor_generalization.R" --cohort AA --run-id "$RUN_ID"

# S-LDSC (null), aging (NOT_SUPPORTED), environmental (exploratory), and the
# schizophrenia locus detail displaced from Figure 5.
log_message "Supplementary figures (Modules 06, 09, 09b, 10) -- AA"
run_r "$HERE/11_supplementary_figures.R" --cohort AA --run-id "$RUN_ID"

# --------------------------------------------------- tables, then the seal
# Table 1 and the cohort QC panels. These used to live in step_2_table1_qc.sh,
# which ran AFTER 03_close_figure_run.R had already set the run read-only --
# so no sealed run ever contained a tables/ directory. They belong here.
log_message "Table 1 and the ancestry panel"
run_r "$HERE/04_table1_cohort.R" --run-id "$RUN_ID"

log_message "Cross-region sample-integrity screen"
run_r "$HERE/05_qc_sample_integrity.R" --run-id "$RUN_ID"

# The registry is built by reading the source_data/ tables the figure builders
# wrote, so it must come after every figure and before the seal.
log_message "Manuscript number registry, claim matrix, denominators, manifest"
run_r "$HERE/10_manuscript_tables.R" --cohort AA --run-id "$RUN_ID"

# Writes the provenance manifest, verifies every figure has source data, and
# seals the run read-only. Must be last.
log_message "Sealing ${RUN_ID}"
run_r "$HERE/03_close_figure_run.R" --run-id "$RUN_ID"

log_message "**** Job ends ****"
