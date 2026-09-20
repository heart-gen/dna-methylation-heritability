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
if (run_is_sealed(manifest)) {
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
    if (is.na(rid) && identical(cell, source_cell)) {
        ## The full-arm run, and ONLY the full-arm run, is keyed without a cell
        ## suffix. Letting a subsample cell take this branch is how the tier
        ## silently compared the arm against itself: the replicate key missed,
        ## the fallback handed back the full caudate run, and every replicate
        ## came out at retained_fraction exactly 1.0 with a NA Wilcoxon p --
        ## which criterion 4 accepts, because 1.0 is finite. A missing
        ## replicate must stop the run, not borrow the number it is supposed to
        ## be tested against.
        rid <- mopt(paste0("upstream_03_local_snp_prediction_", region))
    }
    if (is.na(rid)) {
        stop("Manifest has no Module 03 run for ", cell, " x ", region,
             " (looked for manifest field '", key, "')")
    }
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
    ## `chrom` travels with the locus so the uncertainty block below can be the
    ## chromosome: loci within one are not independent. Fall back to the vmr_id
    ## prefix (chrN:start-end) if the column is ever absent.
    chrom <- if ("chrom" %in% names(dt)) as.character(dt$chrom) else
        sub(":.*$", "", as.character(dt$vmr_id))
    out <- dt[, .(vmr_id = as.character(vmr_id),
                  r2 = suppressWarnings(as.numeric(get(col))))]
    out[, chrom := chrom][, cell := cell][, run_id := rid]
    out[is.finite(r2)]
}

full <- load_r2(source_cell, ds_region)
subs <- rbindlist(lapply(sub_cells, load_r2, region = ds_region),
                  use.names = TRUE)

## The comparison is only meaningful if each replicate is a DIFFERENT run from
## the full arm and from every other replicate. This asserts on the resolved run
## IDs rather than trusting the manifest key shape, so a future key change that
## re-points a replicate at the arm fails here instead of producing a
## retained_fraction of exactly 1.0 that the gate would accept.
resolved <- c(full = full$run_id[[1]],
              vapply(split(subs$run_id, subs$cell), function(x) x[[1]],
                     character(1)))
if (anyDuplicated(resolved)) {
    dup <- resolved[duplicated(resolved) | duplicated(resolved, fromLast = TRUE)]
    stop("Tier 3 resolved the same Module 03 run for more than one cell: ",
         paste(sprintf("%s=%s", names(dup), dup), collapse = ", "),
         ". A donor-count sensitivity cannot compare a run against itself.")
}

## ------------------------------------------- estimator resolution at lower n
##
## PI 2026-09-18: retain the boundary-rate shift as an explicit sample-size
## sensitivity, phrased as an ESTIMATOR-RESOLUTION effect and not as biology.
##
## Module 02's `boundary_rate` is the fraction of eligible loci whose unbounded
## joint-PVE estimate sits at or outside the frozen model's output range. In
## caudate it is entirely the LOWER boundary (zero upper hits in the arm and in
## all three replicates), so it measures the mass of loci with no detectable
## local genetic control, and it rises as donors are removed. The magnitudes are
## NOT quoted here on purpose: they depend on which 02 run the arm resolves to
## (the caudate rescore moved the arm's rate), so they live on the output rows
## and in the run's interpretation-constraints.txt, both derived at run time.
##
## That matters for how this tier is read. Part of any attenuation after
## downsampling is the estimator losing resolution, not the caudate changing.
## The numbers are carried on the output rows so the caveat cannot be separated
## from the result.
load_boundary_rate <- function(cell, region) {
    key <- paste0("upstream_02_local_genetic_variance_",
                  gsub("[.]", "_", cell), "_", region)
    rid <- mopt(key)
    if (is.na(rid) && identical(cell, source_cell)) {
        rid <- mopt(paste0("upstream_02_local_genetic_variance_", region))
    }
    if (is.na(rid)) return(list(run_id = NA_character_, rate = NA_real_))
    f <- file.path(V2_ROOT, "02_local_genetic_variance", "_m", "runs", rid,
                   "results", "combined", "observed-score-decision.tsv")
    if (!file.exists(f)) return(list(run_id = rid, rate = NA_real_))
    d <- as.data.table(fread(f))
    list(run_id = rid,
         rate = if ("boundary_rate" %in% names(d)) {
             suppressWarnings(as.numeric(d$boundary_rate[[1]]))
         } else NA_real_)
}
boundary_full <- load_boundary_rate(source_cell, ds_region)
boundary_sub <- lapply(sub_cells, load_boundary_rate, region = ds_region)
names(boundary_sub) <- sub_cells

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
        ## ---- PRIMARY ENDPOINT (PI 2026-09-18) -------------------------------
        ## The within-caudate change on the identical shared caudate locus set,
        ## paired on locus: full n=153 against this n=118 replicate. Both terms
        ## come from the same region and the same loci, so no cross-region
        ## reference enters and no locus-set asymmetry applies. This -- not the
        ## fraction of the DLPFC gap closed -- is what establishes whether donor
        ## count moves the estimate.
        primary_endpoint = as.character(ds$primary_endpoint),
        mean_delta_r2 = mean(delta),
        ## Positive = the draw-down lowered r2, the direction donor count
        ## predicts. Stated explicitly so the sign is not read backwards.
        attenuation_mean_r2 = mean(a$r2) - mean(b$r2),
        ## A_r = (R2_full - R2_n118,r) / R2_full, the PI-specified effect size.
        ## Identical to 1 - retained_fraction; carried under its own name
        ## because it is the quantity the primary criterion is stated in.
        relative_attenuation = (mean(a$r2) - mean(b$r2)) / mean(a$r2),
        retained_fraction = mean(b$r2) / mean(a$r2),
        spearman_full_vs_subset = stats::cor(a$r2, b$r2, method = "spearman"),
        ## DESCRIPTIVE ONLY. It gates nothing: PI 2026-09-18 replaced the
        ## significance criterion with the effect size A, because with ~11k
        ## paired loci a trivial attenuation would be "significant" and would
        ## still say nothing about whether donor count matters. Kept because it
        ## costs nothing and a reader will look for it.
        paired_wilcoxon_p = wt$p.value,
        paired_wilcoxon_p_is_not_a_gate = TRUE,
        ## ---- estimator resolution, not biology ------------------------------
        boundary_rate_full = boundary_full$rate,
        boundary_rate_subset = boundary_sub[[this_cell]]$rate,
        boundary_rate_delta = boundary_sub[[this_cell]]$rate - boundary_full$rate,
        boundary_shift_is_estimator_resolution_not_biology = TRUE)
}), use.names = TRUE)

write_atomic(per_rep, file.path(out_dir, "caudate-downsampling-replicates.tsv"))

## ------------------------------------- secondary: the region-level gap context
##
## PI 2026-09-18. `closed` is SECONDARY and DESCRIPTIVE. Writing it out:
##
##   closed = 1 - (m_sub - ref)/(m_full - ref) = (m_full - m_sub)/(m_full - ref)
##
## The numerator is the primary within-caudate attenuation -- `ref` cancels, so
## it is clean and locus-matched. The DENOMINATOR is the caudate-minus-DLPFC
## gap, and its DLPFC term is a mean over a DIFFERENT VMR population: caudate and
## DLPFC carry different vmr_set_ids, so no cross-region locus intersection
## exists and the caudate means are restricted to `shared` while the DLPFC mean
## is not. Selection into `shared` is not random with respect to r2, so the
## denominator carries an uncontrolled term.
##
## The ratio therefore contextualizes HOW MUCH OF THE REGION-LEVEL GAP the
## primary attenuation would represent. It does not establish the primary
## result, and it cannot on its own decide whether donor count "explains" the
## excess. The asymmetry is written into interpretation-constraints.txt.
dlpfc_r2 <- tryCatch(load_r2(source_cell, "dlpfc"),
                     error = function(e) NULL)
excess_full <- NA_real_; excess_sub <- NA_real_; closed <- NA_real_
if (!is.null(dlpfc_r2)) {
    ref <- mean(dlpfc_r2$r2)
    excess_full <- mean(full[match(shared, vmr_id)]$r2) - ref
    excess_sub <- mean(per_rep$mean_r2_subset) - ref
    ## Guarded: a near-zero denominator would make this explode and read as a
    ## finding.
    closed <- if (is.finite(excess_full) && abs(excess_full) > 1e-6) {
        1 - excess_sub / excess_full
    } else NA_real_
}

## ------------------------------------------------------------- the reading
##
## Both thresholds now come from the pi_locked config (see
## caudate_downsampling.replicate_agreement_max_range and
## .major_contributor_gap_closed_min). They were literals here until
## 2026-09-18, which left the two numbers that select the reading as the only
## ones in this tier not prespecified -- and editable without touching a locked
## file. Values are unchanged.
agreement_max_range <- as.numeric(ds$replicate_agreement_max_range)
gap_closed_min <- as.numeric(ds$major_contributor_gap_closed_min)
if (!is.finite(agreement_max_range) || !is.finite(gap_closed_min)) {
    stop("caudate_downsampling lacks finite replicate_agreement_max_range / ",
         "major_contributor_gap_closed_min; refusing to fall back to a literal")
}

## ---------------------------------------- uncertainty on A: REPORTED, not gated
##
## PI 2026-09-18: report a chromosome-block CI on the attenuation as uncertainty,
## and do NOT make statistical significance another hard gate. The criteria
## above are an agreement tolerance and two effect-size/interpretive thresholds;
## none of them is a significance cutoff.
##
## Delete-one-chromosome weighted block jackknife (Busing, Meijer & van der
## Leeden 1999), the same construction 06_partitioned_heritability uses for
## block standard errors. Weighted because chromosomes differ several-fold in
## locus count, and the unweighted formula assumes equal blocks.
relative_attenuation_of <- function(keep_ids) {
    a <- full[match(keep_ids, vmr_id)]
    mean(vapply(sub_cells, function(this_cell) {
        b <- subs[cell == this_cell][match(keep_ids, vmr_id)]
        (mean(a$r2) - mean(b$r2)) / mean(a$r2)
    }, numeric(1)))
}
A_hat <- mean(per_rep$relative_attenuation)
shared_chrom <- full[match(shared, vmr_id)]$chrom
blocks <- split(shared, shared_chrom)
blocks <- blocks[lengths(blocks) > 0L]
A_ci_lower <- A_ci_upper <- A_se <- NA_real_
n_blocks <- length(blocks)
if (n_blocks >= 2L) {
    n_tot <- length(shared)
    hj <- n_tot / lengths(blocks)
    theta_j <- vapply(names(blocks), function(b) {
        relative_attenuation_of(setdiff(shared, blocks[[b]]))
    }, numeric(1))
    ## Pseudo-values, then the weighted jackknife variance.
    pseudo <- hj * A_hat - (hj - 1) * theta_j
    theta_J <- mean(pseudo)
    A_se <- sqrt(sum((pseudo - theta_J)^2 / (hj - 1)) / n_blocks)
    z <- stats::qnorm(1 - (1 - as.numeric(ds$attenuation_ci_level)) / 2)
    A_ci_lower <- A_hat - z * A_se
    A_ci_upper <- A_hat + z * A_se
}

## ------------------------------------------------------------- the reading
##
## Three locked criteria, PI 2026-09-18. (i) and (ii) are the PRIMARY and are
## evaluated on the shared caudate loci alone; (iii) is separate and is the only
## one that mentions the region-level difference.
##
##   (i)   all replicates attenuate in the same direction
##   (ii)  mean relative attenuation A >= primary_min_relative_attenuation
##   (iii) gap_closed >= major_contributor_gap_closed_min  -- SEPARATE, and what
##         licenses calling donor count a plausible major contributor to the
##         caudate-DLPFC difference specifically
##
## Plus the agreement tolerance: replicates whose A_r spread exceeds it leave the
## tier unable to resolve the question at this n. Note A_r = 1 - retained_r, so
## the spread is the same number either way.
min_rel_atten <- as.numeric(ds$primary_min_relative_attenuation)
if (!is.finite(min_rel_atten)) {
    stop("caudate_downsampling lacks a finite primary_min_relative_attenuation")
}
agree <- diff(range(per_rep$relative_attenuation)) <= agreement_max_range
## (i) direction
attenuation_consistent <- all(per_rep$relative_attenuation > 0) ||
    all(per_rep$relative_attenuation <= 0)
primary_attenuation <- mean(per_rep$attenuation_mean_r2)
## (ii) magnitude
primary_meets_magnitude <- is.finite(A_hat) && A_hat >= min_rel_atten
primary_established <- attenuation_consistent && primary_meets_magnitude
## (iii) the separate interpretive criterion
gap_criterion_met <- is.finite(closed) && closed >= gap_closed_min

reading <- if (!agree || !attenuation_consistent) {
    "indeterminate_replicates_disagree"
} else if (primary_established && gap_criterion_met) {
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
    ## ---- PRIMARY: within-caudate, shared locus set, no cross-region term ----
    primary_endpoint = as.character(ds$primary_endpoint),
    n_donors_full = per_rep$n_donors_full[[1]],
    mean_r2_full = per_rep$mean_r2_full[[1]],
    mean_r2_subset_mean = mean(per_rep$mean_r2_subset),
    primary_attenuation_mean_r2 = primary_attenuation,
    ## A = mean_r (R2_full - R2_n118,r)/R2_full -- the criterion quantity.
    relative_attenuation_mean = A_hat,
    relative_attenuation_min = min(per_rep$relative_attenuation),
    relative_attenuation_max = max(per_rep$relative_attenuation),
    relative_attenuation_range = diff(range(per_rep$relative_attenuation)),
    ## Uncertainty: REPORTED, never gated (PI 2026-09-18).
    relative_attenuation_se_blockjack = A_se,
    relative_attenuation_ci_lower = A_ci_lower,
    relative_attenuation_ci_upper = A_ci_upper,
    attenuation_uncertainty_method = as.character(ds$attenuation_uncertainty),
    attenuation_ci_level = as.numeric(ds$attenuation_ci_level),
    attenuation_n_blocks = n_blocks,
    attenuation_significance_is_not_a_gate = TRUE,
    ## The three locked criteria, each with the threshold that judged it.
    criterion_i_same_direction = attenuation_consistent,
    criterion_ii_magnitude_met = primary_meets_magnitude,
    primary_min_relative_attenuation = min_rel_atten,
    primary_established = primary_established,
    criterion_iii_gap_closed_met = gap_criterion_met,
    major_contributor_gap_closed_min = gap_closed_min,
    retained_fraction_mean = mean(per_rep$retained_fraction),
    retained_fraction_min = min(per_rep$retained_fraction),
    retained_fraction_max = max(per_rep$retained_fraction),
    replicates_agree = agree,
    ## An agreement TOLERANCE, not a significance cutoff.
    replicate_agreement_max_range = agreement_max_range,
    ## ---- estimator resolution at lower n, NOT biology ----------------------
    boundary_rate_full = boundary_full$rate,
    boundary_rate_subset_mean = mean(per_rep$boundary_rate_subset),
    boundary_rate_delta_mean = mean(per_rep$boundary_rate_delta),
    boundary_shift_is_estimator_resolution_not_biology = TRUE,
    ## Reported ALONGSIDE the attenuation, never folded into it: A is computed on
    ## prediction r2 alone and carries no boundary-rate correction, so overall
    ## prediction attenuation stays distinguishable from the accompanying loss of
    ## estimator resolution (PI 2026-09-18).
    boundary_shift_excluded_from_attenuation_threshold = TRUE,
    ## ---- SECONDARY / DESCRIPTIVE: the region-level gap context -------------
    dlpfc_reference_mean_r2 = if (!is.null(dlpfc_r2)) mean(dlpfc_r2$r2) else NA_real_,
    excess_over_dlpfc_full = excess_full,
    excess_over_dlpfc_subset = excess_sub,
    fraction_of_excess_closed_by_matching_n = closed,
    fraction_of_excess_closed_is_descriptive_only = TRUE,
    ## The denominator's DLPFC term is a mean over a different VMR population,
    ## and the caudate means are restricted to `shared` while it is not.
    gap_ratio_locus_set_asymmetry = paste0(
        "caudate means on ", length(shared), " shared loci (vmr_set_id ",
        "differs from dlpfc); dlpfc reference on ",
        if (!is.null(dlpfc_r2)) nrow(dlpfc_r2) else NA_integer_,
        " dlpfc loci, unrestricted; no cross-region locus intersection exists"),
    reading = reading,
    reading_established_by = "primary within-caudate paired delta, sized by the secondary gap ratio",
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

message("[tier3] PRIMARY -- within-caudate, ", length(shared),
        " shared loci, paired on locus:")
print(per_rep[, .(replicate_cell, mean_r2_full, mean_r2_subset,
                  attenuation_mean_r2, retained_fraction, paired_wilcoxon_p)])
message("[tier3] estimator resolution (not biology):")
print(per_rep[, .(replicate_cell, boundary_rate_full, boundary_rate_subset,
                  boundary_rate_delta)])
message("[tier3] SECONDARY / descriptive -- region-level gap context:")
print(summary_dt[, .(dlpfc_reference_mean_r2, excess_over_dlpfc_full,
                     fraction_of_excess_closed_by_matching_n,
                     fraction_of_excess_closed_is_descriptive_only)])
message("[tier3] locked criteria:")
print(summary_dt[, .(criterion_i_same_direction, relative_attenuation_mean,
                     criterion_ii_magnitude_met, criterion_iii_gap_closed_met,
                     replicates_agree, reading)])
message("[tier3] attenuation A = ", sprintf("%.4f", A_hat),
        if (is.finite(A_se)) sprintf(" (%.0f%% block-jackknife CI %.4f to %.4f, %d chromosome blocks; reported, not gated)",
                                     100 * as.numeric(ds$attenuation_ci_level),
                                     A_ci_lower, A_ci_upper, n_blocks) else "")
message("[tier3] reading: ", reading)
message("  Established by the PRIMARY within-caudate paired comparison on one ",
        "locus set; the DLPFC gap ratio only sizes it and is descriptive.")
message("  Part of any attenuation is ESTIMATOR RESOLUTION at lower n, not a ",
        "biological change: the lower-boundary mass rises from ",
        sprintf("%.4f", boundary_full$rate), " at n=",
        per_rep$n_donors_full[[1]], " to ",
        sprintf("%.4f", mean(per_rep$boundary_rate_subset)), " at n=", target_n, ".")
message("  This tier licenses ONLY whether donor count explains the caudate ",
        "excess. A surviving residual is NOT biological: caudate stays ",
        "confounded with sequencing batch whatever this shows (AGENTS.md 8.1).")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
