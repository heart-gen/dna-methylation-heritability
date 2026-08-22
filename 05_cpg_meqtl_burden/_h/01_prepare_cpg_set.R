#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- assemble the testable CpG set ####
##
## Usage:
##   Rscript _h/01_prepare_cpg_set.R --run-id cmb-AA-caudate-20260817
##
## Produces the CpG x donor methylation matrix, the CpG position file and the
## covariate file that the mapping step consumes, plus -- and this is the part
## the legacy analysis got wrong -- an explicit accounting of every CpG that is
## a member of a VMR but will NOT be tested, with the reason.
##
## AGENTS.md 7.5: "Report tested CpGs separately from prepared-but-untested
## CpGs. Audit every concordance denominator." A burden statistic of the form
## n_cpgs_with_meqtl / n_cpgs is only interpretable if the denominator is the
## TESTED count. Using membership as the denominator silently penalizes VMRs
## whose CpGs failed coverage QC, and coverage correlates with repeat content,
## which is the confound 04 is also fighting.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "05_cpg_meqtl_burden"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mval("cohort"); region <- mval("region")
meqtl <- load_config("meqtl_parameters")

cat_dir <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs",
                     mval("upstream_vmr_catalog_run_id"), "results")

## ------------------------------------------------------- corrected membership
membership <- fread(file.path(cat_dir, "cpg_vmr_membership.tsv"))
if (!all(c("cpg_id", "vmr_id", "chrom", "pos") %in% names(membership))) {
    stop("Membership table is missing required columns; expected ",
         "cpg_id, vmr_id, chrom, pos")
}
if (!identical(as.character(membership$vmr_set_id[1]), mval("vmr_set_id"))) {
    stop("Membership table cites vmr_set_id ", membership$vmr_set_id[1],
         " but this run was opened against ", mval("vmr_set_id"))
}

## ------------------------------------------------------------- methylation
## Non-residualized CpG methylation is the primary source (meqtl_parameters.yml).
## Residualized values are permitted only as a negative-control sensitivity:
## residualizing on PCs derived from the same CpGs removes part of the genetic
## signal being tested.
meth <- fread(file.path(cat_dir, "cpg_meth.phen"))
donors <- as.character(meth[[2]])
assert_no_dups(donors, "donors in the CpG methylation matrix")
assert_expected_n(length(donors), cohort, region)

cpg_ids <- names(meth)[-(1:2)]

## --------------------------------------------------------------- exclusions
## Every exclusion is recorded with a reason, and the reasons are mutually
## exclusive and ordered, so a CpG appears exactly once in the audit.
excl <- data.table(cpg_id = character(), reason = character())
add_excl <- function(ids, reason) {
    ids <- setdiff(ids, excl$cpg_id)
    if (length(ids) > 0) {
        excl <<- rbind(excl, data.table(cpg_id = ids, reason = reason))
    }
}

add_excl(setdiff(membership$cpg_id, cpg_ids), "not_in_methylation_matrix")

qc <- meqtl$cpg_qc
mat <- as.matrix(meth[, -(1:2)])
frac_nonmissing <- colMeans(!is.na(mat))
add_excl(cpg_ids[frac_nonmissing < qc$min_fraction_samples_passing_coverage],
         "below_min_fraction_samples_passing_coverage")

if (isTRUE(qc$require_nonzero_variance)) {
    v <- apply(mat, 2, stats::var, na.rm = TRUE)
    add_excl(cpg_ids[!is.finite(v) | v == 0], "zero_variance")
}

if (isTRUE(qc$exclude_common_ct_snp_at_cpg)) {
    ct_dir <- resolve_path("ct_snp_dir")
    ct_files <- list.files(ct_dir, full.names = TRUE)
    if (length(ct_files) == 0) {
        stop("exclude_common_ct_snp_at_cpg is TRUE but no C>T SNP files under ", ct_dir)
    }
    ct <- unique(rbindlist(lapply(ct_files, fread), fill = TRUE))
    ## Match on position, not on ID: the C>T lists and the BSobj use different
    ## identifier conventions and an ID join silently matches nothing.
    ct_key <- paste(ct[[1]], ct[[2]], sep = ":")
    mem_key <- paste(membership$chrom, membership$pos, sep = ":")
    add_excl(membership$cpg_id[mem_key %in% ct_key], "common_ct_snp_at_cpg")
}

tested <- setdiff(intersect(membership$cpg_id, cpg_ids), excl$cpg_id)
if (length(tested) == 0) stop("No CpGs survived QC; check the exclusion reasons")

## ------------------------------------------------------------------ outputs
out <- file.path(run_dir, "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

write_atomic(membership[cpg_id %in% tested], file.path(out, "tested-cpg-membership.tsv"))
write_atomic(excl, file.path(out, "excluded-cpgs.tsv"))

## The denominator table. Every downstream burden fraction must use
## n_tested_cpgs from here, never n_member_cpgs.
denom <- membership[, .(n_member_cpgs = .N), by = vmr_id]
denom <- merge(denom,
               membership[cpg_id %in% tested, .(n_tested_cpgs = .N), by = vmr_id],
               by = "vmr_id", all.x = TRUE)
denom[is.na(n_tested_cpgs), n_tested_cpgs := 0L]
write_atomic(denom, file.path(out, "vmr-cpg-denominators.tsv"))

summary_dt <- data.table(
    region = region, population = cohort, vmr_set_id = mval("vmr_set_id"),
    n_member_cpgs = nrow(membership),
    n_tested_cpgs = length(tested),
    n_excluded_cpgs = nrow(excl),
    n_vmrs_with_zero_tested = denom[n_tested_cpgs == 0, .N]
)
write_atomic(summary_dt, file.path(out, "cpg-set-summary.tsv"))
write_atomic(excl[, .N, by = reason], file.path(out, "exclusion-reasons.tsv"))
print(summary_dt)

append_manifest(list(dir = run_dir), list(
    n_donors = length(donors),
    donor_checksum = donor_checksum(donors),
    n_tested_cpgs = length(tested)
))
