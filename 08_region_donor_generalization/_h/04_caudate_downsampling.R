#!/usr/bin/env Rscript
#### 08 Stage 04 -- tier 3: does donor count explain the caudate excess? ####
##
## Usage:
##   Rscript _h/04_caudate_downsampling.R --run-id rdg-AA-crossregion-YYYYMMDD
##
## Module 03 attributed caudate's prediction excess to its donor count (153 vs
## DLPFC's 118) and never tested it. This tier tests exactly that and nothing
## else, by comparing the full caudate run against three independent draws of
## 118 of its own donors, refit end to end on the SAME locus set.
##
## Only two readings are permitted, and they are named in the config:
##   * the excess largely disappears -> donor count is a plausible major
##     contributor;
##   * it persists -> donor count does not explain it.
##
## A THIRD reading is explicitly forbidden: a surviving excess does not become
## biological. Caudate remains perfectly confounded with sequencing batch
## whatever this shows (AGENTS.md 8.1), so "not donor count" never becomes
## "therefore region". The output carries that as a column, not as a footnote,
## because this is the tier most likely to be over-read.
##
## Three replicates, not one: a single draw would put the conclusion at the
## mercy of one partition of 153 donors. The replicate spread IS part of the
## result -- if the three draws disagree, the answer is that the design cannot
## resolve it at this n.

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
mopt <- function(f) {
    v <- manifest$value[manifest$field == f]
    if (length(v) != 1L) NA_character_ else as.character(v[[1L]])
}
if (nzchar(manifest$value[manifest$field == "finished_at"][1] %||% "")) {
    stop("Run is sealed and immutable: ", opts$run_id)
}

cfg <- load_config("region_donor_generalization")
out_dir <- file.path(run_dir, "results")
tier <- "mechanistic_sensitivity"
ds <- config_get(cfg, "caudate_downsampling")

if (!isTRUE(ds$enabled)) {
    message("[tier3] caudate_downsampling.enabled is false; nothing to do")
    quit(save = "no", status = 0)
}

ds_region <- as.character(ds$region)
source_cell <- as.character(ds$source_cell)
target_n <- as.integer(ds$target_n)
n_rep <- as.integer(ds$n_replicates)
sub_cells <- paste0(source_cell, ".", sprintf("n%dr%d", target_n, seq_len(n_rep)))

## ------------------------------------------------ the quantity being compared
##
## Module 03's per-locus out-of-fold r2. It is the statistic the attribution was
## made about, and unlike the Module 02 score it is NOT a within-cell rank -- it
## is an r2 on the same phenotype scale in every cell, so its LEVEL is
## comparable across donor counts. That is what makes this tier possible at all.
load_r2 <- function(cell, region) {
    key <- paste0("upstream_03_local_snp_prediction_",
                  gsub("[.]", "_", cell), "_", region)
    rid <- mopt(key)
    if (is.na(rid)) {
        ## The full-arm run is keyed without a cell suffix.
        rid <- mopt(paste0("upstream_03_local_snp_prediction_", region))
    }
    if (is.na(rid)) stop("Manifest has no Module 03 run for ", cell, " x ", region)
    d <- file.path(V2_ROOT, "03_local_snp_prediction", "_m", "runs", rid,
                   "results", "combined")
    f <- list.files(d, pattern = "^oof-prediction.*\\.tsv$", full.names = TRUE)
    f <- f[!grepl("per-donor", f)]
    if (length(f) != 1L) {
        stop("Expected one per-locus OOF table in ", d, ", found ", length(f))
    }
    dt <- as.data.table(fread(f[[1L]]))
    col <- if ("r2_pred_oof" %in% names(dt)) "r2_pred_oof" else
        stop("No r2_pred_oof column in ", f[[1L]])
    if (!"vmr_id" %in% names(dt)) stop("No vmr_id column in ", f[[1L]])
    out <- dt[, .(vmr_id = as.character(vmr_id),
                  r2 = suppressWarnings(as.numeric(get(col))))]
    out[, cell := cell][, run_id := rid]
    out[is.finite(r2)]
}

full <- load_r2(source_cell, ds_region)
subs <- rbindlist(lapply(sub_cells, load_r2, region = ds_region),
                  use.names = TRUE)

## The locus set must be held fixed, or the comparison confounds donor count
## with locus turnover. Restrict to loci scored in the full run AND in every
## replicate, and report what that costs.
shared <- Reduce(intersect, c(list(full$vmr_id),
                              split(subs$vmr_id, subs$cell)))
if (length(shared) < 100L) {
    stop("Only ", length(shared), " loci are scored in the full run and every ",
         "replicate; a donor-count sensitivity needs the locus set held fixed.")
}

## ------------------------------------------------------- per-replicate result
##
## Paired on locus, because the pairing is the whole design: the same locus,
## the same window, the same SNPs available, only the donor count differs.
per_rep <- rbindlist(lapply(sub_cells, function(this_cell) {
    ## Distinct name: `cell` is also a column of `subs`, and inside the
    ## data.table `i` expression the bare name resolves to the COLUMN, which
    ## would match every row. That is defect V2 in miniature.
    a <- full[match(shared, vmr_id)]
    b <- subs[cell == this_cell][match(shared, vmr_id)]
    delta <- b$r2 - a$r2
    wt <- stats::wilcox.test(b$r2, a$r2, paired = TRUE, exact = FALSE)
    data.table(
        tier = tier,
        region = ds_region,
        replicate_cell = this_cell,
        replicate_run = b$run_id[[1]],
        n_loci = length(shared),
        ## From the locked config, never a literal: if the arm's design_n
        ## moves, this must move with it or the tier silently compares the
        ## wrong two numbers.
        n_donors_full = as.integer(
            load_config("cohorts")$donor_counts[[source_cell]][[ds_region]]$design_n),
        n_donors_subset = target_n,
        mean_r2_full = mean(a$r2),
        mean_r2_subset = mean(b$r2),
        median_r2_full = stats::median(a$r2),
        median_r2_subset = stats::median(b$r2),
        pct_positive_full = 100 * mean(a$r2 > 0),
        pct_positive_subset = 100 * mean(b$r2 > 0),
        mean_delta_r2 = mean(delta),
        ## The headline: how much of the full-run mean survives the draw-down.
        retained_fraction = mean(b$r2) / mean(a$r2),
        spearman_full_vs_subset = stats::cor(a$r2, b$r2, method = "spearman"),
        paired_wilcoxon_p = wt$p.value)
}), use.names = TRUE)

write_atomic(per_rep, file.path(out_dir, "caudate-downsampling-replicates.tsv"))

## ------------------------------------------------------------- the reading
##
## The DLPFC reference is what makes "the excess" a defined quantity: the excess
## is caudate-minus-DLPFC, so the question is whether drawing caudate down to
## DLPFC's n closes the gap.
dlpfc_r2 <- tryCatch(load_r2(source_cell, "dlpfc"),
                     error = function(e) NULL)
excess_full <- NA_real_; excess_sub <- NA_real_; closed <- NA_real_
if (!is.null(dlpfc_r2)) {
    ref <- mean(dlpfc_r2$r2)
    excess_full <- mean(full[match(shared, vmr_id)]$r2) - ref
    excess_sub <- mean(per_rep$mean_r2_subset) - ref
    ## Fraction of the original excess removed by matching n. Guarded: a
    ## near-zero denominator would make this explode and read as a finding.
    closed <- if (is.finite(excess_full) && abs(excess_full) > 1e-6) {
        1 - excess_sub / excess_full
    } else NA_real_
}

agree <- diff(range(per_rep$retained_fraction)) <= 0.10
reading <- if (!agree) {
    "indeterminate_replicates_disagree"
} else if (is.finite(closed) && closed >= 0.5) {
    "donor_count_is_a_plausible_major_contributor"
} else {
    "donor_count_does_not_explain_the_excess"
}

permitted <- as.character(ds$permitted_readings)
forbidden <- as.character(ds$forbidden_readings)

## The config states the permitted readings as the PI wrote them, in prose
## ("donor count is a plausible major contributor"); this stage derives a
## snake_case token. Comparing the two vocabularies directly made the check
## below unsatisfiable: BOTH substantive readings failed it and only
## `indeterminate_replicates_disagree` passed, so Stage 04 could complete only
## when the three draws disagreed -- the one outcome that resolves nothing.
## `config/region_donor_generalization.yml` is `pi_locked`, so the prose stays
## exactly as written and is still what lands in the output column; only the
## comparison is normalised.
as_reading_token <- function(x) {
    gsub("_+", "_", gsub("[^a-z0-9]+", "_", tolower(trimws(as.character(x)))))
}
permitted_tokens <- c(as_reading_token(permitted),
                      "indeterminate_replicates_disagree")

summary_dt <- data.table(
    tier = tier,
    region = ds_region,
    licenses = trimws(config_get(cfg, paste0("tiers.", tier, ".licenses"))),
    n_replicates = n_rep,
    n_loci_held_fixed = length(shared),
    n_loci_full_run = nrow(full),
    target_n = target_n,
    mean_r2_full = per_rep$mean_r2_full[[1]],
    mean_r2_subset_mean = mean(per_rep$mean_r2_subset),
    retained_fraction_mean = mean(per_rep$retained_fraction),
    retained_fraction_min = min(per_rep$retained_fraction),
    retained_fraction_max = max(per_rep$retained_fraction),
    replicates_agree = agree,
    dlpfc_reference_mean_r2 = if (!is.null(dlpfc_r2)) mean(dlpfc_r2$r2) else NA_real_,
    excess_over_dlpfc_full = excess_full,
    excess_over_dlpfc_subset = excess_sub,
    fraction_of_excess_closed_by_matching_n = closed,
    reading = reading,
    ## Carried on the row so the constraint cannot be separated from the number.
    permitted_readings = paste(permitted, collapse = "; "),
    forbidden_readings = paste(forbidden, collapse = "; "),
    residual_excess_may_be_called_biological = FALSE,
    caudate_remains_batch_confounded = TRUE)
write_atomic(summary_dt, file.path(out_dir, "caudate-downsampling-summary.tsv"))

if (!as_reading_token(reading) %in% permitted_tokens) {
    stop("Stage 04 produced the reading '", reading,
         "', which is not in caudate_downsampling.permitted_readings (",
         paste(permitted, collapse = "; "), ")")
}

print(per_rep[, .(replicate_cell, mean_r2_full, mean_r2_subset,
                  retained_fraction, paired_wilcoxon_p)])
print(summary_dt[, .(retained_fraction_mean, replicates_agree,
                     fraction_of_excess_closed_by_matching_n, reading)])
message("[tier3] reading: ", reading)
message("  This tier licenses ONLY whether donor count explains the caudate ",
        "excess. A surviving residual is NOT biological: caudate stays ",
        "confounded with sequencing batch whatever this shows (AGENTS.md 8.1).")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
