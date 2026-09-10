#!/usr/bin/env Rscript
#### 02 / stage 14: recombine the donor-group estimation cells ####
##
## Joins the two donor-group cells of one region -- all_individuals.AA and
## all_individuals.EA -- into one per-locus table with _AA / _EA suffixed
## columns. This is the v1 recombination
## (local-snp-prediction/.../venn_diagram/_h/02.compare_cohort.R:53-62, an inner
## join on locus coordinates with group suffixes) rebuilt on the v2 endpoint.
##
## WHAT THIS DELIBERATELY DOES NOT DO
##
## It emits no pooled rank, no cross-cell score difference and no combined PVE.
## `local_snp_contribution_score` is a within-cell midrank percentile;
## config/local_genetic_control.yml locks `rank_scope: cohort_by_region` and
## AGENTS.md 7.6 forbids raw score-level comparison across cells. A pooled rank
## over two cells would be exactly the prohibited quantity, and a difference of
## two percentiles computed in different denominators is not an effect size.
##
## What the table IS for: asking whether the ORDERING of loci agrees between
## donor groups (rank correlation, concordance of quartile membership), with the
## per-cell n, num_snps, p_eff and ld_metric sitting beside every row so that
## sample size, MAF, LD and SNP availability can be eliminated before any
## difference is discussed (AGENTS.md 7.7).
##
## Usage:
##   Rscript _h/14_combine_donor_group_cells.R --region dlpfc \
##     --run-ids lgv-all_individuals.AA-dlpfc-20260910,lgv-all_individuals.EA-dlpfc-20260910

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

opts <- parse_v2_args(require = c("region", "run_ids"))
region <- tolower(opts$region)
validate_cohort_region(region = region)

run_ids <- trimws(strsplit(opts$run_ids, ",", fixed = TRUE)[[1]])
run_ids <- run_ids[nzchar(run_ids)]
if (length(run_ids) != 2L) {
    stop("Give exactly two run IDs, one per donor group: --run-ids a,b")
}

module_root <- file.path(V2_ROOT, "02_local_genetic_variance")

## Columns that are per-cell measurements and get suffixed. Everything else is a
## property of the LOCUS and must be identical between cells; the join asserts
## that rather than silently keeping one side's copy.
CELL_COLS <- c(
    "n", "num_snps", "n_variants", "snps_in_window", "p_eff", "ld_metric",
    "mean_methylation", "methylation_variance",
    "bslmm_pve", "he_h2", "he_se", "he_pvalue", "rho2_oof", "r2_oof",
    "feature_complete", "computational_failure", "terminal_status",
    "exclusion_reason", "joint_pve_domain_status",
    "pve_cis_joint_unbounded", "pve_cis_joint_calibrated",
    "local_genetic_control_eligible", "local_genetic_control_exclusion_reason",
    "local_snp_contribution_score", "local_snp_contribution_score_z",
    "local_snp_contribution_quartile"
)
LOCUS_COLS <- c("vmr_id", "chrom", "start", "end", "n_cpgs", "vmr_set_id")

read_cell <- function(run_id) {
    run_dir <- file.path(module_root, "_m", "runs", run_id)
    manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
    mval <- function(f) {
        v <- manifest$value[manifest$field == f]
        if (length(v) != 1L) stop("Run ", run_id, " manifest lacks ", f)
        as.character(v[[1L]])
    }
    cell <- mval("cohort")
    parsed <- parse_cell(cell)
    if (!parsed$is_estimation_cell) {
        stop("Run ", run_id, " is discovery arm '", cell, "'. A donor-group ",
             "contrast uses two estimation cells on ONE pooled catalog; ",
             "contrasting an arm against its own superset is the nested ",
             "comparison AGENTS.md 7.7 forbids.")
    }
    if (!identical(tolower(mval("region")), region)) {
        stop("Run ", run_id, " is region '", mval("region"), "', not ", region)
    }
    if (!identical(toupper(mval("smoke_run")), "FALSE")) {
        stop("Run ", run_id, " is a smoke run and is not citable")
    }
    ## The gate, not a convenience: only an ACCEPTED cell may be recombined.
    accepted <- read_accepted_runs("02_local_genetic_variance")
    if (!nrow(accepted) || !run_id %in% accepted$run_id) {
        stop("Run ", run_id, " is not recorded in the 02 README's accepted ",
             "runs table (AGENTS.md 6)")
    }

    tbl <- load_local_genetic_control(run_id, region = region, cohort = cell,
                                      eligible_only = FALSE)
    tbl <- as.data.table(tbl)
    missing <- setdiff(c(LOCUS_COLS, CELL_COLS), names(tbl))
    if (length(missing)) {
        stop("Run ", run_id, " table lacks: ", paste(missing, collapse = ", "))
    }
    list(run_id = run_id, cell = cell, group = parsed$estimation_group,
         catalog_cohort = parsed$catalog_cohort,
         vmr_set_id = mval("vmr_set_id"),
         upstream_catalog = mval("upstream_catalog_run_id"),
         tbl = tbl[, c(LOCUS_COLS, CELL_COLS), with = FALSE])
}

cells <- lapply(run_ids, read_cell)
names(cells) <- vapply(cells, function(x) x$group, character(1))

if (anyDuplicated(names(cells))) {
    stop("Both runs estimate in the same donor group (",
         names(cells)[[1]], "). That is not a contrast.")
}
if (length(unique(vapply(cells, function(x) x$catalog_cohort, character(1)))) != 1L) {
    stop("The two cells were discovered on different catalogs; their loci are ",
         "not the same set.")
}
## The whole design in one assertion: discovery happened ONCE.
if (length(unique(vapply(cells, function(x) x$vmr_set_id, character(1)))) != 1L) {
    stop("The two cells carry different vmr_set_ids. Discovery must happen ",
         "once in the pooled sample (AGENTS.md 7.7); refusing to join.")
}
if (length(unique(vapply(cells, function(x) x$upstream_catalog, character(1)))) != 1L) {
    stop("The two cells trace to different Module 01 catalog runs")
}

## ------------------------------------------------------------------- join
a <- cells[[1]]; b <- cells[[2]]
locus_a <- a$tbl[, LOCUS_COLS, with = FALSE][order(vmr_id)]
locus_b <- b$tbl[, LOCUS_COLS, with = FALSE][order(vmr_id)]
if (!identical(nrow(locus_a), nrow(locus_b)) || !isTRUE(all.equal(locus_a, locus_b))) {
    stop("The two cells' locus tables differ. On a shared catalog every VMR ",
         "must appear in both cells with identical coordinates.")
}

suffix_cell <- function(x) {
    t <- copy(x$tbl)
    setnames(t, CELL_COLS, paste0(CELL_COLS, "_", x$group))
    t[, (setdiff(LOCUS_COLS, "vmr_id")) := NULL]
    t
}

joined <- merge(locus_a, suffix_cell(a), by = "vmr_id", all = FALSE)
joined <- merge(joined, suffix_cell(b), by = "vmr_id", all = FALSE)
if (nrow(joined) != nrow(locus_a)) {
    stop("Join lost ", nrow(locus_a) - nrow(joined), " loci; the two cells ",
         "should agree on every VMR")
}

ga <- a$group; gb <- b$group
elig_a <- joined[[paste0("local_genetic_control_eligible_", ga)]]
elig_b <- joined[[paste0("local_genetic_control_eligible_", gb)]]
elig_a <- as.logical(elig_a); elig_b <- as.logical(elig_b)

## Loci eligible in only one cell are KEPT, with the reason recorded. Dropping
## them would silently condition the comparison on joint eligibility, which is
## itself a function of the very things (n, MAF, SNP availability) the donor
## groups differ in.
joined[, comparable_in_both := elig_a & elig_b]
joined[, comparability := fifelse(
    elig_a & elig_b, "both",
    fifelse(elig_a, paste0("only_", ga),
            fifelse(elig_b, paste0("only_", gb), "neither")))]

joined[, `:=`(
    region = region,
    catalog_cohort = a$catalog_cohort,
    vmr_set_id = a$vmr_set_id,
    donor_groups = paste(ga, gb, sep = ","),
    upstream_lgv_run_id_1 = a$run_id,
    upstream_lgv_run_id_2 = b$run_id,
    upstream_catalog_run_id = a$upstream_catalog,
    ## Carried on every row, as Module 02's own tables do. The rank basis is
    ## unchanged by joining two cells, and neither cell's score became absolute.
    absolute_pve_interpretation_allowed = FALSE,
    rank_scope = "cohort_by_region",
    pooled_rank_emitted = FALSE
)]

out_dir <- file.path(module_root, "_m", "combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(out_dir,
                   paste0("local-genetic-control-donor-group-", region, ".tsv"))
write_atomic(joined, out_f)

## Ordering agreement, the quantity this table exists to support. Spearman on
## the within-cell percentiles is invariant to the fact that the two ranks were
## computed in different denominators, which is exactly why it is reportable
## where a raw difference is not.
both <- joined[comparable_in_both == TRUE]
score_a <- both[[paste0("local_snp_contribution_score_", ga)]]
score_b <- both[[paste0("local_snp_contribution_score_", gb)]]
summary_dt <- data.table(
    region = region,
    n_loci_shared = nrow(joined),
    n_loci_comparable = nrow(both),
    n_only_a = sum(joined$comparability == paste0("only_", ga)),
    n_only_b = sum(joined$comparability == paste0("only_", gb)),
    n_neither = sum(joined$comparability == "neither"),
    group_a = ga, group_b = gb,
    n_donors_a = both[[paste0("n_", ga)]][1],
    n_donors_b = both[[paste0("n_", gb)]][1],
    median_num_snps_a = stats::median(both[[paste0("num_snps_", ga)]], na.rm = TRUE),
    median_num_snps_b = stats::median(both[[paste0("num_snps_", gb)]], na.rm = TRUE),
    median_p_eff_a = stats::median(both[[paste0("p_eff_", ga)]], na.rm = TRUE),
    median_p_eff_b = stats::median(both[[paste0("p_eff_", gb)]], na.rm = TRUE),
    spearman_score = if (nrow(both) > 2)
        stats::cor(score_a, score_b, method = "spearman", use = "complete.obs")
        else NA_real_,
    quartile_concordance = if (nrow(both) > 0) mean(
        both[[paste0("local_snp_contribution_quartile_", ga)]] ==
        both[[paste0("local_snp_contribution_quartile_", gb)]], na.rm = TRUE)
        else NA_real_,
    absolute_pve_interpretation_allowed = FALSE
)
write_atomic(summary_dt, file.path(
    out_dir, paste0("donor-group-concordance-", region, ".tsv")))

print(summary_dt)
message("[recombine] ", nrow(joined), " shared loci -> ", out_f)
message("  Ordering agreement only. Do NOT read the two cells' scores as ",
        "comparable levels, and do not attribute any difference to ancestry ",
        "before eliminating n, MAF, LD and SNP availability (AGENTS.md 7.7).")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
