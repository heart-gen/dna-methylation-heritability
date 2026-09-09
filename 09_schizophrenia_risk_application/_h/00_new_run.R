#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA --region caudate [--allow-unlocked]
##
## Six upstreams are consumed: 01 for the corrected VMR boundaries, 02 for the
## relative local-control score, 04 for the repeat/repressive features the
## integration analysis tests against, 05 for the nominal CpG meQTL statistics
## that link risk variants to VMRs, 06 for the S-LDSC context that decides
## whether this module adds anything, and 07 for transcriptional coupling. All
## must describe the same vmr_set_id, or the module would join disease
## evidence, architecture and coupling for different loci sharing an ID.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "09_schizophrenia_risk_application"
MODULE_TAG <- "scz"

opts <- parse_v2_args(require = c("cohort", "region"))
allow_unlocked <- isTRUE(opts$allow_unlocked)

scz <- load_config("schizophrenia")
assert_locked(list(schizophrenia = scz), allow_unlocked = allow_unlocked)

## AGENTS.md 7.8 and 3: the disease application is tested against the relative
## local SNP contribution score. Legacy predictability (r_squared_cv, defect E1)
## and absolute PVE are banned, and the ban is checked before a directory
## exists rather than after a week of compute.
if (!identical(scz$testing$architecture_predictor,
               "local_snp_contribution_score_z")) {
    stop("config/schizophrenia.yml must use local_snp_contribution_score_z as ",
         "the architecture predictor; legacy predictability and absolute h2 ",
         "are banned (AGENTS.md 3).")
}
if (!isTRUE(scz$locus_definition$independent_of_methylation)) {
    stop("locus_definition.independent_of_methylation must be true: PGC loci ",
         "are defined without reference to any methylation result ",
         "(AGENTS.md 7.8).")
}
if (!isTRUE(scz$testing$separate_fdr_family)) {
    stop("testing.separate_fdr_family must be true: risk-variant x CpG tests ",
         "keep their own FDR family (AGENTS.md 7.8).")
}
if (is.null(scz$prioritization$rule) ||
    !identical(scz$prioritization$rule$name, "ranked_composite_v1")) {
    stop("prioritization.rule is unset or unrecognised. The rule must be ",
         "prespecified before results are seen (AGENTS.md 7.8).")
}

## Prespecified external inputs must exist before a run directory is created.
for (key in c("pgc3_loci_hg38", "pgc3_index_snp_xls", "pgc3_sumstats_hg19",
              "liftover_chain_hg19_to_hg38")) {
    f <- scz$locus_definition[[key]]
    if (is.null(f) || !file.exists(f)) {
        stop("Missing prespecified input '", key, "': ", f %||% "<null>")
    }
}

## Only arms the config enables are opened, and the enabled set is frozen into
## the manifest here so an arm cannot be added after results are seen.
arms <- names(scz$colocalization$arms)[
    vapply(scz$colocalization$arms, function(a) isTRUE(a$enabled), logical(1))]
if (length(arms) == 0) stop("No colocalization arm is enabled in config")
gate_arms <- names(scz$colocalization$arms)[
    vapply(scz$colocalization$arms, function(a) isTRUE(a$gate_eligible), logical(1))]
if (length(gate_arms) == 0) {
    stop("No colocalization arm is gate-eligible. At least one ancestry-matched ",
         "arm is required for the coloc gate to mean anything (AGENTS.md 7.8).")
}

upstreams <- list(
    vmr_catalog             = "01_vmr_catalog",
    local_genetic_variance  = "02_local_genetic_variance",
    repeat_architecture     = "04_repeat_repressive_architecture",
    cpg_meqtl_burden        = "05_cpg_meqtl_burden",
    partitioned_h2          = "06_partitioned_heritability",
    tx_coupling             = "07_transcription_splicing_coupling"
)
accepted <- lapply(upstreams, require_accepted_upstream,
                   cohort = opts$cohort, region = opts$region,
                   allow_unaccepted = allow_unlocked)

sets <- unlist(lapply(accepted, function(a) a$vmr_set_id))
sets <- unique(sets[!is.na(sets)])
if (length(sets) > 1) {
    stop("vmr_set_id mismatch across upstreams: ", paste(sets, collapse = " vs "),
         "\n  01, 02, 04, 05, 06 and 07 must describe the same VMR set ",
         "(AGENTS.md 6).")
}

if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked); ",
         "production run IDs are derived, not chosen.")
}

run <- new_run(
    module = MODULE_TAG, cohort = opts$cohort, region = opts$region,
    module_root = file.path(repo_root(), MODULE),
    run_id = opts$run_id,
    vmr_set_id = if (length(sets)) sets[1] else NA_character_,
    upstream = stats::setNames(
        lapply(accepted, function(a) a$run_id %||% NA_character_),
        paste0(names(accepted), "_run_id")),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_schizophrenia_sha256 = attr(scz, "config_sha256"),
        primary_region = scz$primary_region,
        fdr_family = scz$testing$fdr_family,
        fdr_alpha = scz$testing$fdr_alpha,
        coloc_arms = paste(arms, collapse = ","),
        coloc_gate_arms = paste(gate_arms, collapse = ","),
        coloc_method = scz$colocalization$method,
        prioritization_rule = scz$prioritization$rule$name,
        ## Module 08 owns the caudate downsampling arm and does not exist. This
        ## field is what 12_apply_gates.R reads to record
        ## caudate_not_sample_size_artifact as PENDING_MODULE_08 rather than
        ## silently passing or failing it.
        module_08_downsampling_available =
            if (isTRUE(scz$gates$require_module_08_downsampling)) "TRUE" else "FALSE"
    )
)

for (d in c("results", "results/figures", "coloc", "coloc/regions", "gwas")) {
    dir.create(file.path(run$dir, d), showWarnings = FALSE, recursive = TRUE)
}
message("[09] run ", run$run_id, " opened (coloc arms: ",
        paste(arms, collapse = ", "), "; gate-eligible: ",
        paste(gate_arms, collapse = ", "), ")")
cat(run$run_id, "\n", sep = "")
