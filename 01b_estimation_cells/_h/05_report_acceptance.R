#!/usr/bin/env Rscript
#### 01b / stage 05: independently re-check a sealed cell's acceptance gate ####
##
## Stage 04 enforces these criteria before it will seal a run. This stage checks
## them AGAIN, from the sealed artifacts rather than from the in-flight state,
## and prints the evidence a human needs to record the run in README.md.
##
## Re-checking is not redundant. Stage 04 asserts as it goes and then makes the
## directory read-only; this reads what was actually written, so a criterion
## that passed on a value held in memory but was recorded wrong on disk is
## caught. It also emits the markdown row itself, so the acceptance table cannot
## drift from the evidence behind it.
##
## This stage NEVER edits README.md. Recording an accepted run is a human act
## (AGENTS.md 6); this only supplies the evidence.
##
## Usage:
##   Rscript _h/05_report_acceptance.R                     # every sealed run
##   Rscript _h/05_report_acceptance.R --run-id <id>       # just one

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

opts <- if (length(commandArgs(trailingOnly = TRUE))) {
    parse_v2_args(require = character())
} else list()

module_root <- file.path(V2_ROOT, "01b_estimation_cells")
runs_root <- file.path(module_root, "_m", "runs")
run_ids <- if (!is.null(opts$run_id)) opts$run_id else
    sort(list.files(runs_root, pattern = "^estcell-"))
if (!length(run_ids)) stop("No 01b runs found under ", runs_root)

cohorts_cfg <- load_config("cohorts")

check_run <- function(run_id) {
    run_dir <- file.path(runs_root, run_id)
    fail <- character(0)
    note <- function(...) fail <<- c(fail, paste0(...))

    if (!file.exists(file.path(run_dir, "output_checksums.tsv"))) {
        return(list(run_id = run_id, sealed = FALSE, pass = FALSE,
                    fail = "not sealed"))
    }
    manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
    mv <- function(f) {
        v <- manifest$value[manifest$field == f]
        if (length(v) != 1L) NA_character_ else as.character(v[[1L]])
    }
    cell <- mv("cohort"); region <- mv("region")
    group <- mv("estimation_group"); catalog <- mv("catalog_cohort")
    source_run <- mv("upstream_vmr_catalog")
    source_vmr <- file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs", source_run,
                            "vmr")

    ## -- 1. every VMR accounted for, zero unaccounted, zero failures
    n_vmrs <- as.integer(mv("n_vmrs"))
    n_ext <- as.integer(mv("n_vmrs_extracted"))
    n_nosnp <- as.integer(mv("n_vmrs_no_cis_variant"))
    if (is.na(n_vmrs) || is.na(n_ext) || is.na(n_nosnp)) {
        note("manifest lacks extraction counts")
    } else if (n_ext + n_nosnp != n_vmrs) {
        note("VMRs unaccounted: ", n_vmrs - n_ext - n_nosnp)
    }
    rec_f <- file.path(run_dir, "task_reconciliation.tsv")
    if (!file.exists(rec_f)) note("no task_reconciliation.tsv") else {
        rec <- fread(rec_f)
        get <- function(k) {
            v <- rec$n[rec$category == k]
            if (length(v) == 1L) as.integer(v) else NA_integer_
        }
        for (k in c("failed", "unaccounted", "unexpected")) {
            if (!identical(get(k), 0L)) note(k, " = ", get(k))
        }
        if (!identical(get("expected"), n_vmrs)) note("reconciliation expected != n_vmrs")
    }

    ## -- 2. locus set checksum-matches the source catalog
    for (f in c("vmr.bed", "vmr_catalog.tsv")) {
        a <- file.path(source_vmr, f); b <- file.path(run_dir, "vmr", f)
        if (!file.exists(a) || !file.exists(b)) note("missing ", f) else if
            (!identical(file_sha256(a), file_sha256(b))) note(f, " differs from catalog")
    }

    ## -- 3. strict subset; cells on this catalog disjoint and exhaustive
    donors <- fread(file.path(run_dir, "vmr", "donors_plink.txt"),
                    header = FALSE, colClasses = "character")[[1]]
    pooled <- fread(file.path(source_vmr, "donors_plink.txt"),
                    header = FALSE, colClasses = "character")[[1]]
    if (!all(donors %in% pooled)) note("donors are not a subset of the pooled set")
    if (length(donors) >= length(pooled)) note("not a strict subset")
    siblings <- Filter(function(nm) {
        identical(cohorts_cfg$estimation_cells[[nm]]$catalog_cohort, catalog)
    }, names(cohorts_cfg$estimation_cells))
    ## Build the sibling run IDs from parts. sub() does not vectorise over
    ## `replacement`, so swapping the cell token with a vector silently used
    ## only the first sibling.
    stamp <- sub(paste0("^estcell-", cell, "-", region, "-"), "", run_id)
    sib_dirs <- file.path(runs_root,
                          paste0("estcell-", siblings, "-", region, "-", stamp))
    names(sib_dirs) <- siblings
    sib_donors <- lapply(sib_dirs, function(d) {
        f <- file.path(d, "vmr", "donors_plink.txt")
        if (file.exists(f)) fread(f, header = FALSE, colClasses = "character")[[1]]
        else NULL
    })
    if (all(vapply(sib_donors, Negate(is.null), logical(1)))) {
        flat <- unlist(sib_donors, use.names = FALSE)
        if (anyDuplicated(flat)) note("donor groups overlap")
        if (!setequal(flat, pooled)) note("donor groups do not cover the pooled set")
    } else {
        note("sibling cell not built; disjointness unverified")
    }

    ## -- 4. observed n matches the config reconstruction
    expected_n <- cohorts_cfg$donor_counts[[cell]][[region]]$phenotype_n
    if (is.null(expected_n)) note("no donor_counts entry") else if
        (length(donors) != expected_n) {
            note("n = ", length(donors), ", config says ", expected_n)
        }

    ## -- 5. PCs cover exactly the donors, in order, non-degenerate
    pc_f <- file.path(run_dir, "covs", "genotype_pcs.tsv")
    n_pc <- NA_integer_
    if (!file.exists(pc_f)) note("no genotype_pcs.tsv") else {
        pcs <- fread(pc_f, colClasses = list(character = c("FID", "IID")))
        pc_cols <- grep("^snpPC[0-9]+$", names(pcs), value = TRUE)
        n_pc <- length(pc_cols)
        if (!identical(as.character(pcs$FID), donors)) note("PC donors differ or are misordered")
        if (!n_pc) note("no snpPC columns") else {
            for (nm in pc_cols) {
                v <- as.numeric(pcs[[nm]])
                if (anyNA(v)) note(nm, " has missing values")
                else if (stats::sd(v) == 0) note(nm, " has zero variance")
            }
        }
        wanted <- length(load_config("covariates")$estimation_cells$genotype_pcs)
        if (n_pc != wanted) note("expected ", wanted, " PCs, found ", n_pc)
    }

    list(run_id = run_id, sealed = TRUE, pass = length(fail) == 0L,
         fail = if (length(fail)) paste(fail, collapse = "; ") else "",
         cell = cell, region = region, group = group,
         vmr_set_id = mv("vmr_set_id"), source_run = source_run,
         n_donors = length(donors), n_vmrs = n_vmrs,
         n_extracted = n_ext, n_nosnp = n_nosnp, n_pc = n_pc)
}

results <- lapply(run_ids, check_run)

cat("\n=== 01b acceptance evidence ===\n\n")
for (r in results) {
    cat(sprintf("%-46s %s\n", r$run_id,
                if (!r$sealed) "NOT SEALED"
                else if (r$pass) "ALL FIVE PASS" else paste("FAIL:", r$fail)))
    if (r$sealed) {
        cat(sprintf("    %s/%s  donors=%d  VMRs=%d (extracted %d, no_cis %d)  PCs=%d\n",
                    r$cell, r$region, r$n_donors, r$n_vmrs, r$n_extracted,
                    r$n_nosnp, r$n_pc))
    }
}

pass <- Filter(function(r) isTRUE(r$pass), results)
if (length(pass)) {
    cat("\n=== README rows for the runs that pass ===\n\n")
    for (r in pass) {
        cat(sprintf("| %s | %s | %s | %s | %s |  | ALL_FIVE_CRITERIA_PASS | %d donors, %d VMRs (%d no cis variant), %d within-group PCs; loci from %s |\n",
                    r$run_id, r$cell, r$region, r$vmr_set_id,
                    format(Sys.Date()), r$n_donors, r$n_vmrs, r$n_nosnp,
                    r$n_pc, r$source_run))
    }
    cat("\nThe accepted_by column is deliberately blank: recording an accepted\n",
        "run is a human act (AGENTS.md 6).\n", sep = "")
}
invisible(NULL)
