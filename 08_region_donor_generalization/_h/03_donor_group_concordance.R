#!/usr/bin/env Rscript
#### 08 Stage 03 -- the donor-group axis: concordance, never difference ####
##
## Usage:
##   Rscript _h/03_donor_group_concordance.R --run-id rdg-AA-crossregion-YYYYMMDD
##
## Reads the already-computed recombination outputs and assembles the axis's one
## reportable statement. It does not re-rank, re-fit or pool anything: the
## Module 02 score is a WITHIN-CELL midrank percentile
## (config/local_genetic_control.yml locks rank_scope: cohort_by_region), so
## there is no pooled quantity to compute and AGENTS.md 7.6 forbids inventing
## one.
##
## The reportable quantity is ORDERING AGREEMENT on the shared locus set. Three
## things make that a real statement rather than a hedge:
##
##   * discovery happened ONCE, in the pooled sample, so neither group's VMR
##     calling can advantage it (AGENTS.md 7.7);
##   * the two donor sets are disjoint and partition the pooled arm exactly;
##   * agreement is read against the ANALYTIC RELIABILITY CEILING from
##     02 Stage 16. Without the ceiling an imperfect rho reads as a donor-group
##     difference when most of it is input uncertainty -- the single most likely
##     misreading of this axis.
##
## Loci eligible in only one cell are retained and labelled. Dropping them would
## condition the comparison on joint eligibility, which is itself a function of
## n, MAF and SNP availability -- the very things the cells differ in.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "08_region_donor_generalization"

opts <- parse_v2_args(require = c("run_id"))
run_dir <- file.path(V2_ROOT, MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]
    if (length(v) != 1L) stop("Run manifest lacks unique field: ", f)
    as.character(v[[1L]])
}
if (run_is_sealed(manifest)) {
    stop("Run is sealed and immutable: ", opts$run_id)
}

cfg <- load_config("region_donor_generalization")
out_dir <- file.path(run_dir, "results")

## The policy is re-read here, not inherited from the manifest. If someone
## widened the config after the run was opened, this stage must refuse rather
## than produce a table under a policy the run never gated on.
policy <- donor_group_inference_policy()
if (!identical(policy$mode, mval("donor_group_inference"))) {
    stop("The donor-group policy changed after this run was opened: manifest ",
         "records '", mval("donor_group_inference"), "', config now says '",
         policy$mode, "'. Open a new run.")
}

cells <- strsplit(mval("donor_group_cells"), ",", fixed = TRUE)[[1]]
groups <- vapply(cells, function(c) parse_cell(c)$estimation_group, character(1))
dg_regions <- as.character(config_get(cfg, "donor_group.regions"))
inputs <- config_get(cfg, "donor_group.inputs")
use_ceiling <- isTRUE(config_get(cfg,
                                 "donor_group.interpret_against_reliability_ceiling"))

resolve <- function(key, region) {
    file.path(V2_ROOT, gsub("{region}", region, inputs[[key]], fixed = TRUE))
}

ga <- groups[[1]]; gb <- groups[[2]]

## ------------------------------------------------------------- per region
rows <- rbindlist(lapply(dg_regions, function(re) {
    lgc <- as.data.table(fread(resolve("local_genetic_control", re)))
    con <- as.data.table(fread(resolve("concordance_summary", re)))

    ## Refuse a stale recombination output. The upstream run IDs in our manifest
    ## are the ones this run gated on; if the combined table was rebuilt from
    ## different runs since, its numbers do not belong to this run.
    for (i in 1:2) {
        f <- paste0("upstream_lgv_run_id_", i)
        if (f %in% names(lgc)) {
            cited <- unique(as.character(lgc[[f]]))
            if (length(cited) != 1L) {
                stop("Combined table for ", re, " cites ", length(cited),
                     " values of ", f)
            }
            key <- paste0("upstream_02_local_genetic_variance_",
                          gsub("[.]", "_", cells[[i]]), "_", re)
            expected <- manifest$value[manifest$field == key]
            if (length(expected) == 1L && !identical(cited, expected)) {
                stop("Combined table for ", re, " cites ", cited,
                     " but this run gated on ", expected,
                     ". Re-run 02/_h/14_combine_donor_group_cells.R.")
            }
        }
    }

    ## Guards that the recombination stage already applies, re-asserted because
    ## this is the stage whose output becomes a manuscript sentence.
    if ("pooled_rank_emitted" %in% names(lgc) &&
        any(toupper(as.character(lgc$pooled_rank_emitted)) == "TRUE")) {
        stop("Combined table for ", re, " emits a pooled rank (AGENTS.md 7.6)")
    }
    if ("absolute_pve_interpretation_allowed" %in% names(lgc) &&
        any(toupper(as.character(lgc$absolute_pve_interpretation_allowed)) == "TRUE")) {
        stop("Combined table for ", re, " authorizes absolute-PVE reading")
    }

    sa <- paste0("local_snp_contribution_score_", ga)
    sb <- paste0("local_snp_contribution_score_", gb)
    qa <- paste0("local_snp_contribution_quartile_", ga)
    qb <- paste0("local_snp_contribution_quartile_", gb)
    for (col in c(sa, sb, qa, qb, "comparable_in_both", "comparability")) {
        if (!col %in% names(lgc)) {
            stop("Combined table for ", re, " lacks ", col)
        }
    }
    ## The column-selection trap that produced a spurious rho = 1.0000 once:
    ## `_z` columns are a monotone transform of their own score, so grepping for
    ## "score" matches four columns and two of them are the same variable.
    ## Name them exactly and assert they are not the same vector.
    both <- lgc[toupper(as.character(comparable_in_both)) == "TRUE"]
    x <- as.numeric(both[[sa]]); y <- as.numeric(both[[sb]])
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    if (length(x) < 3L) stop("Fewer than 3 comparable loci in ", re)
    if (identical(x, y)) {
        stop("The two cells' score columns are identical in ", re,
             ". That is a column-selection error, not a result.")
    }

    rho <- stats::cor(x, y, method = "spearman")
    ## Top-quartile overlap: of the loci in the top quartile in one cell, the
    ## fraction also in the top quartile in the other. Symmetric by taking the
    ## Jaccard, so the statistic does not depend on which cell is named first.
    top_a <- both[[qa]] == max(both[[qa]], na.rm = TRUE)
    top_b <- both[[qb]] == max(both[[qb]], na.rm = TRUE)
    top_overlap <- sum(top_a & top_b, na.rm = TRUE) /
        max(1L, sum(top_a | top_b, na.rm = TRUE))

    out <- data.table(
        tier = "donor_group_concordance",
        region = re,
        group_a = ga, group_b = gb,
        cell_a = cells[[1]], cell_b = cells[[2]],
        n_loci_shared = nrow(lgc),
        n_loci_comparable = length(x),
        n_only_a = lgc[comparability == paste0("only_", ga), .N],
        n_only_b = lgc[comparability == paste0("only_", gb), .N],
        n_neither = lgc[comparability == "neither", .N],
        spearman_score = rho,
        quartile_concordance = as.numeric(con$quartile_concordance[[1]]),
        top_quartile_overlap = top_overlap,
        n_donors_a = as.numeric(con$n_donors_a[[1]]),
        n_donors_b = as.numeric(con$n_donors_b[[1]]))

    ## ------------------------------------------ against the ceiling
    if (use_ceiling) {
        cf <- resolve("reliability_ceiling", re)
        if (!file.exists(cf)) {
            stop("interpret_against_reliability_ceiling is true but the ",
                 "ceiling for ", re, " is missing: ", cf,
                 "\n  Run 02/_h/16_reliability_ceiling.R.")
        }
        ceil <- as.data.table(fread(cf))
        need <- c("ceiling_bslmm", "observed_cross_cell_spearman",
                  "concordance_fraction_of_ceiling_bslmm")
        if (!all(need %in% names(ceil))) {
            stop("Reliability-ceiling table for ", re, " lacks: ",
                 paste(setdiff(need, names(ceil)), collapse = ", "))
        }
        out[, `:=`(
            reliability_ceiling = as.numeric(ceil$ceiling_bslmm[[1]]),
            ceiling_basis = as.character(ceil$ceiling_basis[[1]]),
            fraction_of_ceiling =
                as.numeric(ceil$concordance_fraction_of_ceiling_bslmm[[1]]))]
        ## The ceiling is an upper bound, so observed agreement above it means
        ## the bound was computed from different inputs than the observation.
        if (is.finite(out$reliability_ceiling) &&
            rho > out$reliability_ceiling + 1e-6) {
            stop("Observed Spearman (", signif(rho, 4), ") exceeds the ",
                 "analytic ceiling (", signif(out$reliability_ceiling, 4),
                 ") in ", re, ". The ceiling and the observation are not from ",
                 "the same runs.")
        }
    }
    out
}), use.names = TRUE, fill = TRUE)

## The policy travels with the result.
rows[, `:=`(
    donor_group_inference = policy$mode,
    cross_group_raw_score_comparison = policy$cross_group_raw_score_comparison,
    ancestry_effect_claim_allowed = policy$ancestry_effect_claim_allowed,
    pooled_rank_emitted = FALSE,
    pooled_r2_emitted = FALSE,
    batch_confounded = isTRUE(config_get(cfg, "donor_group.batch_confounded")),
    required_eliminations = paste(
        config_get(cfg, "interpretation.required_eliminations_before_ancestry_language"),
        collapse = ","))]

## Donor counts must reconstruct the pooled arm exactly, or the two cells are
## not a partition and "disjoint sets on a shared locus set" is false.
cohorts <- load_config("cohorts")
for (re in dg_regions) {
    pooled <- cohorts$donor_counts[[config_get(cfg, "donor_group.catalog_cohort")]][[re]]$design_n
    got <- rows[region == re, n_donors_a + n_donors_b]
    if (!identical(as.integer(got), as.integer(pooled))) {
        stop("Donor counts for ", re, " sum to ", got, ", not the locked ",
             pooled, ". The cells do not partition the pooled arm.")
    }
}

write_atomic(rows, file.path(out_dir, "donor-group-concordance.tsv"))

print(rows[, .(region, n_loci_comparable, spearman_score, top_quartile_overlap,
               reliability_ceiling, fraction_of_ceiling)])
message("[donor-group] concordance on the shared pooled-discovery locus set, ",
        "read against the analytic ceiling.")
message("  Reportable: ORDERING AGREEMENT. Not reportable: any level, any ",
        "cross-cell difference, any ancestry attribution ",
        "(config/analysis_thresholds.yml:donor_group; AGENTS.md 7.6, 7.7).")
message("  The cells differ in n, MAF spectrum, LD and SNP availability; all ",
        "must be eliminated before a difference is discussed.")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
