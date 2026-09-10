#!/usr/bin/env Rscript
#### 10_environmental_exploratory -- build and gate the exposure matrix ####
##
## Usage:
##   Rscript _h/01_build_exposure_matrix.R --run-id env-AA-caudate-YYYYMMDD
##
## Two jobs, and the second is the scientific one.
##
## Build: assemble the donor-level exposure table on the donor set the accepted
## Module 01 run actually used, deriving the composites declared in
## config/environmental.yml.
##
## Gate: decide which exposures are testable in THIS cohort x region, before any
## association is computed. This is the prespecified form of the power
## limitation that demoted the analysis from primary to exploratory. v1 tested
## cocaine on two positive donors in caudate and one in dlpfc; a "hit" there is
## a statement about one person's methylome.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "10_environmental_exploratory", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "10_environmental_exploratory"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mf("cohort")
region <- mf("region")
catalog_cohort <- if (is.na(mf("catalog_cohort"))) cohort else mf("catalog_cohort")
env <- load_run_config("environmental", run_dir)

## ---------------------------------------------------------------------------
## Donor set: taken from the accepted Module 01 run, never re-derived.
## ---------------------------------------------------------------------------
## Re-deriving it from the phenotype table would silently disagree with the run
## that defined the VMRs the moment a filter drifted. The legacy Table 1 is
## wrong for exactly this reason (11_integrated_manuscript_outputs/README.md).
vmr_run_id <- mf("upstream_vmr_catalog_run_id")
if (is.na(vmr_run_id)) stop("Manifest carries no upstream_vmr_catalog_run_id")
vmr_run_dir <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs", vmr_run_id)
donor_f <- file.path(vmr_run_dir, "vmr", "donors_plink.txt")
if (!file.exists(donor_f)) stop("Missing Module 01 donor list: ", donor_f)
donors <- fread(donor_f, header = FALSE, colClasses = "character")
setnames(donors, c("brnum", "iid")[seq_len(ncol(donors))])
assert_no_dups(donors$brnum, "Module 01 donors")
assert_expected_n(nrow(donors), cohort, region)

## ---------------------------------------------------------------------------
## Phenotype table
## ---------------------------------------------------------------------------
## cohort_def() already resolves phenotype_table against the repo root.
arm <- cohort_def(catalog_cohort)
pheno_f <- arm$phenotype_table
if (!file.exists(pheno_f)) stop("Missing phenotype table: ", pheno_f)
pheno <- fread(pheno_f, na.strings = c("NA", ""))

needed <- unique(c("brnum", "region", "agedeath", "race",
                   unlist(lapply(env$exposures$binary, `[[`, "source_columns")),
                   as.character(env$exposures$descriptive),
                   vapply(env$exposures$categorical, `[[`, character(1),
                          "source_column")))
missing_cols <- setdiff(needed, names(pheno))
if (length(missing_cols)) {
    stop("Phenotype table lacks declared exposure column(s): ",
         paste(missing_cols, collapse = ", "))
}

## `..region` is column-selection syntax and does not resolve in `i`; a local
## name that differs from the column name is the reliable form.
region_filter <- region
pheno <- pheno[region == region_filter]
pheno <- unique(pheno, by = "brnum")
setkey(pheno, brnum)
missing_donors <- setdiff(donors$brnum, pheno$brnum)
if (length(missing_donors)) {
    stop("Donors in the accepted Module 01 run are absent from the phenotype ",
         "table: ", paste(head(missing_donors, 10), collapse = ", "))
}
## Order follows the Module 01 donor list, so the ordered-donor checksum in the
## manifest describes the same sequence the upstream run recorded (AGENTS.md 9).
dt <- pheno[donors$brnum]

## ---------------------------------------------------------------------------
## Composites
## ---------------------------------------------------------------------------
## One missingness convention everywhere: all source columns NA -> NA, otherwise
## the logical OR. `any(na.rm = TRUE)` over an all-NA row returns FALSE, which
## would convert "not recorded" into "not exposed" and inflate every denominator.
as_logical_col <- function(x) {
    if (is.logical(x)) return(x)
    v <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(v))
    out[v %in% c("true", "t", "1", "yes")] <- TRUE
    out[v %in% c("false", "f", "0", "no")] <- FALSE
    if (any(!is.na(v) & is.na(out))) {
        stop("Unparseable logical value(s): ",
             paste(unique(v[!is.na(v) & is.na(out)]), collapse = ", "))
    }
    out
}
or_composite <- function(dt, cols) {
    m <- vapply(cols, function(cc) as_logical_col(dt[[cc]]), logical(nrow(dt)))
    m <- matrix(m, nrow = nrow(dt))
    all_na <- rowSums(!is.na(m)) == 0
    out <- rowSums(m, na.rm = TRUE) > 0
    out[all_na] <- NA
    out
}

out <- data.table(brnum = dt$brnum, primarydx = as.character(dt$primarydx))
for (nm in names(env$exposures$binary)) {
    spec <- env$exposures$binary[[nm]]
    out[[nm]] <- or_composite(dt, unlist(spec$source_columns))
}
for (nm in as.character(env$exposures$descriptive)) {
    out[[nm]] <- as_logical_col(dt[[nm]])
}
for (nm in names(env$exposures$categorical)) {
    spec <- env$exposures$categorical[[nm]]
    src <- as.character(dt[[spec$source_column]])
    lvl <- rep(NA_character_, length(src))
    for (target in names(spec$collapse)) {
        lvl[src %in% unlist(spec$collapse[[target]])] <- target
    }
    ## A value present in the data but in no collapse bucket becomes NA
    ## silently, which is how a whole education stratum disappears. v1's map
    ## omitted "MD" for exactly this reason.
    unmapped <- setdiff(unique(src[!is.na(src)]),
                        unlist(spec$collapse, use.names = FALSE))
    if (length(unmapped)) {
        stop("config/environmental.yml collapse map for '", nm,
             "' does not cover observed value(s): ",
             paste(unmapped, collapse = ", "))
    }
    out[[nm]] <- factor(lvl, levels = names(spec$collapse))
}

## ---------------------------------------------------------------------------
## Eligibility gate
## ---------------------------------------------------------------------------
n_donors <- nrow(out)

## Thresholds differ by stratum. The pooled floor was locked before any count
## was seen; the stratified floor was lowered to 15 by the PI on 2026-09-10
## after the stratified counts were computed, and config/environmental.yml says
## so in as many words. Read them from config rather than branching on a
## literal, so the provenance stays in one auditable place.
gate_for <- function(stratum) {
    if (identical(stratum, "all")) {
        list(min_minor = as.integer(env$eligibility$min_minor_class_n),
             max_miss = as.numeric(env$eligibility$max_missing_frac))
    } else {
        list(min_minor = as.integer(env$eligibility$stratified$min_minor_class_n),
             max_miss = as.numeric(env$eligibility$stratified$max_missing_frac))
    }
}

## An exposure may declare the strata it is testable in. antipsychotics does:
## there are zero exposed controls, so a pooled test would present a within-case
## contrast as a population-level exposure effect.
declared_strata <- function(nm) {
    spec <- env$exposures$binary[[nm]]
    if (!is.null(spec) && !is.null(spec$strata)) as.character(unlist(spec$strata))
    else names(env$strata)
}

assess <- function(nm, kind, stratum) {
    g <- gate_for(stratum)
    min_minor <- g$min_minor; max_miss <- g$max_miss
    keep_rows <- if (identical(stratum, "all")) rep(TRUE, n_donors) else
        out$primarydx %in% as.character(unlist(env$strata[[stratum]]$dx_filter))
    sub <- out[keep_rows]
    n_stratum <- nrow(sub)
    v <- sub[[nm]]
    nonmiss <- sum(!is.na(v))
    miss_frac <- if (n_stratum == 0) 1 else 1 - nonmiss / n_stratum
    counts <- table(v[!is.na(v)])
    minor <- if (length(counts) < 2) 0L else as.integer(min(counts))
    reasons <- character(0)
    if (length(counts) < 2) reasons <- c(reasons, "zero_variance")
    if (minor < min_minor) {
        reasons <- c(reasons, sprintf("minor_class_%d_below_%d", minor, min_minor))
    }
    if (miss_frac > max_miss) {
        reasons <- c(reasons, sprintf("missing_%.3f_above_%.3f", miss_frac, max_miss))
    }
    if (!stratum %in% declared_strata(nm)) {
        reasons <- c(reasons, "not_declared_for_this_stratum")
    }
    data.table(
        exposure = nm, stratum = stratum, kind = kind,
        n_donors = n_stratum, gate_min_minor_class_n = min_minor,
        gate_max_missing_frac = max_miss,
        n_nonmissing = nonmiss, missing_frac = round(miss_frac, 4),
        minor_class_n = minor,
        class_counts = paste(sprintf("%s=%d", names(counts), as.integer(counts)),
                             collapse = ", "),
        eligible = kind != "descriptive" && length(reasons) == 0,
        reason = if (length(reasons)) paste(reasons, collapse = ";") else
            if (kind == "descriptive") "descriptive_only_union_is_tobacco" else "",
        ## Honesty about the boundary. Both criteria are round numbers, and
        ## several exposures land just the wrong side of one of them:
        ## any_trauma_hx misses on 16.2% missing against a 15% ceiling. A
        ## variable that close is reported as marginal, not as a clean
        ## exclusion, so a reader can see the decision was a threshold rather
        ## than an absence of signal.
        near_gate = kind != "descriptive" && length(reasons) > 0 &&
            !("zero_variance" %in% reasons) &&
            !("not_declared_for_this_stratum" %in% reasons) &&
            minor >= (min_minor - 3L) && miss_frac <= (max_miss + 0.03)
    )
}

elig <- rbindlist(unlist(lapply(names(env$strata), function(st) c(
    lapply(names(env$exposures$binary), assess, kind = "binary", stratum = st),
    lapply(as.character(env$exposures$descriptive), assess,
           kind = "descriptive", stratum = st),
    lapply(names(env$exposures$categorical), assess,
           kind = "categorical", stratum = st))), recursive = FALSE))
elig[, `:=`(cohort = cohort, region = region, run_id = opts$run_id)]

## Variables retired in config are not silently absent from the record.
retired <- rbindlist(lapply(names(env$exposures$retired), function(nm) {
    spec <- env$exposures$retired[[nm]]
    data.table(exposure = nm, stratum = NA_character_,
               kind = "retired_in_config", n_donors = n_donors,
               gate_min_minor_class_n = NA_integer_,
               gate_max_missing_frac = NA_real_,
               n_nonmissing = NA_integer_, missing_frac = NA_real_,
               minor_class_n = NA_integer_,
               class_counts = as.character(spec$counts %||% ""),
               eligible = FALSE, reason = as.character(spec$reason),
               near_gate = FALSE, cohort = cohort, region = region,
               run_id = opts$run_id)
}))

write_atomic(rbind(elig, retired),
             file.path(run_dir, "results", "exposure-eligibility.tsv"))
write_atomic(out, file.path(run_dir, "results", "exposure-matrix.tsv"))

pairs <- elig[eligible == TRUE, .(exposure, stratum)]
if (nrow(pairs) < as.integer(env$gates$min_eligible_exposures)) {
    stop("No exposure clears the eligibility gate in any stratum of ", cohort,
         " x ", region, ". Nothing to test; see results/exposure-eligibility.tsv.")
}
## The unit downstream is the exposure x stratum PAIR: each is its own BH
## family, and the same exposure can be eligible in one stratum and not another.
eligible <- pairs[, paste(exposure, stratum, sep = "@")]

append_manifest(list(dir = run_dir), list(
    n_donors_exposure_matrix = as.character(n_donors),
    donor_checksum_exposure_matrix = donor_checksum(out$brnum),
    eligible_exposures = paste(eligible, collapse = ","),
    n_eligible_exposures = as.character(length(eligible)),
    strata = paste(names(env$strata), collapse = ","),
    ineligible_exposures = paste(
        elig[eligible == FALSE & kind != "descriptive",
             paste(exposure, stratum, sep = "@")], collapse = ","),
    near_gate_exposures = paste(
        elig[near_gate == TRUE, paste(exposure, stratum, sep = "@")],
        collapse = ",")
))

message("[10] ", cohort, " x ", region, ": ", n_donors, " donors; eligible ",
        "exposure@stratum: ", paste(eligible, collapse = ", "))
if (any(elig$near_gate)) {
    message("[10] within 3 donors of the gate (report as marginal, not clean): ",
            paste(elig[near_gate == TRUE,
                       paste(exposure, stratum, sep = "@")], collapse = ", "))
}
