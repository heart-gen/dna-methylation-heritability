#!/usr/bin/env Rscript
#### 08 Stage 06 -- seal the run ####
##
## Usage:
##   Rscript _h/06_finalize_run.R --run-id rdg-AA-crossregion-YYYYMMDD
##
## Checksums every output, records the decision in the manifest and makes the
## directory read-only. Sealing is NOT acceptance: it says the computation
## finished and its outputs are pinned. Acceptance is a PI act recorded in the
## module README (AGENTS.md 6), and 02 Stage 16 is the cautionary case -- it
## once required sealed and not accepted, so a superseded run's numbers looked
## citable.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "08_region_donor_generalization"

opts <- parse_v2_args(require = c("run_id"))
run_dir <- file.path(V2_ROOT, MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
if (nzchar(manifest$value[manifest$field == "finished_at"][1] %||% "")) {
    stop("Run is already sealed: ", opts$run_id)
}

decision_f <- file.path(run_dir, "results",
                        "region-donor-generalization-decision.tsv")
if (!file.exists(decision_f)) {
    stop("No decision table; run _h/05_apply_gates.R before sealing")
}
decision <- as.data.table(fread(decision_f))

## Every stage's output must be present. Sealing a run missing a stage would
## produce a directory that looks complete and is not.
required <- c("tiers.tsv", "cross-region-tests.tsv",
              "cross-region-replication.tsv",
              "cross-region-rank-agreement.tsv", "cross-region-summary.tsv",
              "identified-difference.tsv", "identified-difference-summary.tsv",
              "descriptive-confounded-regions.tsv",
              "donor-group-concordance.tsv",
              "region-donor-generalization-qc.tsv",
              "region-donor-generalization-decision.tsv")
tier3 <- identical(
    toupper(as.character(manifest$value[manifest$field == "tier3_enabled"][1])),
    "TRUE")
if (tier3) {
    required <- c(required, "caudate-downsampling-replicates.tsv",
                  "caudate-downsampling-summary.tsv")
}
missing <- required[!file.exists(file.path(run_dir, "results", required))]
if (length(missing)) {
    stop("Refusing to seal an incomplete run; missing: ",
         paste(missing, collapse = ", "))
}

## The interpretation constraints are written into the run directory as prose,
## the way 04, 05 and 07 do, so a reader who opens only the run gets them.
writeLines(c(
    paste0("08_region_donor_generalization -- ", opts$run_id),
    "",
    "Every result in this run carries exactly one tier. The tier is the only",
    "interpretation the result supports. results/tiers.tsv is the table.",
    "",
    "TIER 1, cross-region replication -- the primary deliverable. Read it as",
    "robustness of the genetic-control architecture across technical AND",
    "regional contexts. Brain region is perfectly confounded with sequencing",
    "batch (AGENTS.md 8.1), so a successful replication is more compelling for",
    "crossing that boundary. Do NOT write that the confounding strengthens a",
    "result: it does nothing for the interpretability of a difference.",
    "",
    "TIER 2, DLPFC vs hippocampus -- the only identified region difference.",
    "Both regions sit inside batches 1-2, so this contrast does not cross the",
    "batch boundary. Differences are claimed only under strict conjunction.",
    "",
    "TIER 3, caudate at n=118 -- licenses whether donor count explains the",
    "caudate excess, and NOTHING more. A surviving residual is not biological.",
    "Caudate remains batch-confounded whatever the downsampling shows, so",
    "'not donor count' never becomes 'therefore region'.",
    "",
    "  TIER 3 PRIMARY ENDPOINT (PI 2026-09-18). The primary result is the",
    "  WITHIN-CAUDATE change on the identical shared caudate locus set: full",
    "  n=153 against each n=118 replicate, paired on locus. Both terms are the",
    "  same region and the same loci, so no cross-region reference enters it.",
    "",
    "  THE GAP RATIO IS SECONDARY AND DESCRIPTIVE.",
    "  fraction_of_excess_closed_by_matching_n equals",
    "      (mean_r2_full - mean_r2_subset) / (mean_r2_full - dlpfc_reference).",
    "  The numerator is the primary within-caudate attenuation and is clean.",
    "  The DENOMINATOR is not: caudate and DLPFC carry DIFFERENT vmr_set_ids, so",
    "  no cross-region locus intersection exists. The caudate means are",
    "  restricted to the loci shared by the full run and all replicates; the",
    "  DLPFC mean is over all DLPFC-scored loci, unrestricted. Selection into",
    "  the shared set is not random with respect to r2, so the denominator",
    "  carries an uncontrolled term. Use the ratio ONLY to say how much of the",
    "  region-level gap the primary attenuation would represent. Do not use it",
    "  on its own to decide whether donor count explains the excess.",
    "",
    "  ESTIMATOR RESOLUTION, NOT BIOLOGY. Module 02's boundary_rate -- the",
    "  fraction of eligible loci whose unbounded estimate sits at the frozen",
    "  model's output floor, i.e. loci with no detectable local genetic control",
    "  -- rises with the draw-down: 0.6263 at n=153 to 0.643-0.646 across the",
    "  three n=118 replicates, spread 0.0025. In caudate it is entirely the",
    "  LOWER boundary; there are zero upper-boundary hits in the arm or any",
    "  replicate. Removing 35 donors therefore pushes about 1.8% more loci below",
    "  the floor. ANY ATTENUATION AFTER DOWNSAMPLING PARTLY REFLECTS STATISTICAL",
    "  RESOLUTION RATHER THAN A BIOLOGICAL CHANGE, and must be reported that",
    "  way. The per-replicate numbers are on the tier-3 output rows.",
    "",
    "TIER 4, caudate vs other regions -- descriptive only. Fitted, written and",
    "surfaced in dedicated columns; supports no claim. Retained because readers",
    "will ask.",
    "",
    "DONOR-GROUP AXIS -- concordance only. The reportable quantity is ORDERING",
    "AGREEMENT between all_individuals.AA and all_individuals.EA on the shared",
    "pooled-discovery locus set, read against the analytic reliability ceiling.",
    "No pooled rank, no pooled r2, no cross-cell score difference, and no",
    "ancestry attribution: the cells differ in n, MAF spectrum, LD and SNP",
    "availability, and AGENTS.md 7.7 requires all of those be eliminated first.",
    "Do NOT contrast AA against all_individuals -- those are nested.",
    "",
    "Donor-group language: AA = 'Black American', EA = 'non-Hispanic white",
    "American'.",
    "",
    paste0("Gate decision: ", decision$decision[[1]], " (",
           decision$n_passed[[1]], "/", decision$n_criteria[[1]],
           " criteria)."),
    "This gate certifies tiering and interpretation constraints, not a positive",
    "finding. Sealing is not acceptance; nothing here is citable until the PI",
    "records this run in the module README (AGENTS.md 6)."
), file.path(run_dir, "results", "interpretation-constraints.txt"))

run <- list(run_id = opts$run_id, dir = run_dir,
            module_root = file.path(V2_ROOT, MODULE))
append_manifest(run, list(
    decision = decision$decision[[1]],
    n_criteria_passed = decision$n_passed[[1]],
    n_criteria = decision$n_criteria[[1]]))
sums <- close_run(run)

message("[seal] ", opts$run_id, ": ", nrow(sums), " outputs, decision ",
        decision$decision[[1]])
message("  Not accepted. Record the run in ", MODULE,
        "/README.md under '## Accepted runs' after PI review (AGENTS.md 6).")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
