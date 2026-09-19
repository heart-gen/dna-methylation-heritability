#!/usr/bin/env Rscript
#### 09 stage 16 -- does the axis contrast survive matching caudate to n=118? ####
##
## Usage:
##   Rscript _h/16_axis_downsampling_sensitivity.R --cohort AA
##   Rscript _h/16_axis_downsampling_sensitivity.R --cohort AA --allow-unlocked
##
## PI 2026-09-18. Module 08 tier 3 showed that donor count is a plausible major
## contributor to caudate's larger PREDICTION magnitude. It did not ask what
## happens to the SCZ axis contrast, which is what decision 2 rests on. Without
## this stage, decision 2 assumes the answer.
##
## Read carefully what this does and does not establish:
##
##   * Tier 3 attenuating says the MAGNITUDE of caudate's local-genetic-control
##     signal depends on donor count. It does NOT say the caudate estimate is
##     biased, and this stage must never be written up that way.
##   * This stage asks a different question: with caudate cut to 118 donors, do
##     SCZ-linked VMRs still sit LOWER on the axis? If the direction holds, the
##     axis contrast is not an artifact of caudate's larger n, and decision 2
##     stands on its own evidence rather than on assumption.
##
## The comparison is legitimate because the axis test is WITHIN cell: Module 02's
## score is a within-cell midrank percentile, so comparing SCZ-linked against
## background inside one cell is exactly what the score supports, and no level is
## ever compared across cells (AGENTS.md 7.6).
##
## SCZ linkage is POSITIONAL -- it comes from published PGC3 intervals and the
## corrected VMR coordinates, neither of which depends on which donors were used
## to fit the SNP model. So the caudate arm's linkage column is reused verbatim
## and only the score is swapped. That is the whole design: one locus set, one
## linkage assignment, three donor draws.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09_schizophrenia_risk_application", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "09_schizophrenia_risk_application"
opts <- parse_v2_args(require = "cohort")
allow_unaccepted <- isTRUE(opts$allow_unlocked)

scz <- load_config("schizophrenia")
rdg <- load_config("region_donor_generalization")
ds <- config_get(rdg, "caudate_downsampling")

predictor <- scz$testing$architecture_predictor
covs_want <- unlist(scz$testing$architecture_covariates)
alpha <- as.numeric(scz$testing$fdr_alpha)
fdr_method <- scz$testing$fdr_method

## The cells come from the SAME locked config Module 08 tier 3 used, so the two
## analyses cannot drift onto different replicate sets.
region <- as.character(ds$region)
source_cell <- as.character(ds$source_cell)
target_n <- as.integer(ds$target_n)
n_rep <- as.integer(ds$n_replicates)
cells <- paste0(source_cell, ".", sprintf("n%dr%d", target_n, seq_len(n_rep)))

module_root <- file.path(repo_root(), MODULE)
out_dir <- file.path(module_root, "_m", "combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------- the full-arm caudate run
full_run <- if (!allow_unaccepted) {
    require_accepted_upstream(MODULE, opts$cohort, region)$run_id
} else {
    acc <- tryCatch(require_accepted_upstream(MODULE, opts$cohort, region)$run_id,
                    error = function(e) NA_character_)
    if (!is.na(acc)) acc else {
        cand <- list.files(file.path(module_root, "_m", "runs"),
                           pattern = paste0("^scz-", opts$cohort, "-", region, "-"))
        if (length(cand) == 0) stop("No Module 09 run for ", region)
        sort(cand, decreasing = TRUE)[[1L]]
    }
}
full_dir <- file.path(module_root, "_m", "runs", full_run)
link_f <- file.path(full_dir, "results", "vmr-scz-linkage.tsv")
if (!file.exists(link_f)) stop("Missing linkage table: ", link_f)
linkage <- as.data.table(fread(link_f))
for (need in c("vmr_id", "scz_linked", predictor)) {
    if (!need %in% names(linkage)) stop("Linkage table lacks '", need, "': ", link_f)
}
linkage[, scz_linked := as.logical(scz_linked)]

man <- fread(file.path(full_dir, "manifest.tsv"), colClasses = "character")
mf <- function(f) {
    v <- man$value[man$field == f]
    if (length(v) != 1L) NA_character_ else as.character(v[[1L]])
}

## Module 04's feature matrix -- joined, not recomputed, so GC content and
## mappability are byte-identical to the published values, exactly as stage 04
## does it.
feat_f <- file.path(repo_root(), "04_repeat_repressive_architecture", "_m",
                    "runs", mf("upstream_repeat_architecture_run_id"),
                    "results", "vmr-features.tsv")
if (!file.exists(feat_f)) stop("Missing Module 04 features: ", feat_f)
feat <- as.data.table(fread(feat_f))
banned <- intersect(unlist(scz$forbidden_columns), names(feat))
if (length(banned)) {
    stop("Module 04 feature table carries superseded estimator column(s): ",
         paste(banned, collapse = ", "))
}
covs <- intersect(covs_want, names(feat))

## -------------------------------------------------------- the axis test itself
##
## Same two models, same covariates, same direction convention as
## _h/04_architecture_axis.R. Kept as one function so the full arm and every
## replicate go through identical code -- a separate code path for the
## sensitivity is how a sensitivity ends up not comparable to its primary.
axis_test <- function(score_dt, label, cell, run_id, n_donors) {
    dt <- merge(linkage[, .(vmr_id, scz_linked)], score_dt, by = "vmr_id")
    dt <- merge(dt, feat[, c("vmr_id", covs), with = FALSE], by = "vmr_id",
                all.x = TRUE)
    score <- dt[[predictor]]
    ok <- is.finite(score)
    if (sum(ok & dt$scz_linked) < 2L || sum(ok & !dt$scz_linked) < 2L) {
        stop("Too few scored VMRs in ", label, " to test the axis contrast")
    }
    w <- suppressWarnings(stats::wilcox.test(score[ok & dt$scz_linked],
                                             score[ok & !dt$scz_linked],
                                             alternative = "two.sided"))
    fit_dt <- dt[ok & complete.cases(dt[, c(predictor, covs), with = FALSE])]
    est <- se <- pv <- NA_real_
    if (nrow(fit_dt) > 0 && length(unique(fit_dt$scz_linked)) == 2L) {
        form <- stats::as.formula(paste("scz_linked ~", predictor,
            if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""))
        co <- summary(stats::glm(form, data = fit_dt,
                                 family = stats::binomial()))$coefficients
        if (predictor %in% rownames(co)) {
            est <- co[predictor, "Estimate"]
            se  <- co[predictor, "Std. Error"]
            pv  <- co[predictor, "Pr(>|z|)"]
        }
    }
    data.table(
        comparison = label, cell = cell, region = region,
        local_genetic_variance_run = run_id, n_donors = n_donors,
        n_vmrs_scored = sum(ok),
        n_linked = sum(ok & dt$scz_linked),
        n_background = sum(ok & !dt$scz_linked),
        unadj_mean_diff = mean(score[ok & dt$scz_linked]) -
                          mean(score[ok & !dt$scz_linked]),
        unadj_p = w$p.value,
        adj_log_odds_per_sd = est, adj_se = se, adj_p = pv,
        adj_or_per_sd = exp(est),
        adj_or_lower = exp(est - 1.96 * se),
        adj_or_upper = exp(est + 1.96 * se),
        covariates = paste(covs, collapse = ","))
}

score_of <- function(cell) {
    lgv <- require_accepted_upstream("02_local_genetic_variance", cell, region,
                                    allow_unaccepted = allow_unaccepted)
    d <- load_local_genetic_control(lgv$run_id, region = region, cohort = cell)
    if (!predictor %in% names(d)) {
        stop("Module 02 table for ", cell, " has no '", predictor, "' column")
    }
    man_c <- fread(file.path(repo_root(), "02_local_genetic_variance", "_m",
                             "runs", lgv$run_id, "manifest.tsv"),
                   colClasses = "character")
    nd <- man_c$value[man_c$field == "n_donors"]
    list(dt = d[, c("vmr_id", predictor), with = FALSE],
         run_id = lgv$run_id,
         n_donors = if (length(nd) == 1L) as.integer(nd) else NA_integer_)
}

full_score <- score_of(source_cell)
rows <- list(axis_test(full_score$dt, "full_arm", source_cell,
                       full_score$run_id, full_score$n_donors))
for (ce in cells) {
    s <- score_of(ce)
    rows[[length(rows) + 1L]] <- axis_test(s$dt, "n_matched_replicate", ce,
                                           s$run_id, s$n_donors)
}
res <- rbindlist(rows, use.names = TRUE)

## One FDR family over the replicate tests, applied once. The full arm is the
## reference, not a member of the family.
res[, fdr_family := "scz_axis_caudate_downsampling"]
res[comparison == "n_matched_replicate",
    adj_q := p.adjust(adj_p, method = fdr_method)]
res[comparison == "full_arm", adj_q := NA_real_]
res[, significant := is.finite(adj_q) & adj_q < alpha]
res[, direction := fifelse(!is.finite(adj_log_odds_per_sd), NA_character_,
                    fifelse(adj_log_odds_per_sd > 0, "higher_in_scz_linked",
                            "lower_in_scz_linked"))]

full_row <- res[comparison == "full_arm"]
reps <- res[comparison == "n_matched_replicate"]

direction_preserved <- all(reps$direction == full_row$direction[[1]], na.rm = TRUE)
n_sig <- sum(reps$significant %in% TRUE)
## The magnitude change, reported on the SAME footing as Module 08 tier 3's A so
## the two are directly comparable: a positive value means the contrast shrank.
rel_change <- (abs(full_row$adj_log_odds_per_sd[[1]]) -
               abs(reps$adj_log_odds_per_sd)) /
              abs(full_row$adj_log_odds_per_sd[[1]])

verdict <- if (!direction_preserved) {
    "axis_contrast_direction_does_not_survive_n_matching"
} else if (n_sig == nrow(reps)) {
    "axis_contrast_survives_n_matching_in_every_replicate"
} else if (n_sig >= 1L) {
    "axis_contrast_direction_survives_significance_partially"
} else {
    "axis_contrast_direction_survives_but_loses_significance"
}

summ <- data.table(
    cohort = opts$cohort, region = region, endpoint = "scz_axis_under_n_matching",
    predictor = predictor,
    n_donors_full = full_row$n_donors[[1]], n_donors_matched = target_n,
    n_replicates = nrow(reps),
    full_or_per_sd = full_row$adj_or_per_sd[[1]],
    full_direction = full_row$direction[[1]],
    replicate_or_per_sd_mean = mean(reps$adj_or_per_sd),
    replicate_or_per_sd_min = min(reps$adj_or_per_sd),
    replicate_or_per_sd_max = max(reps$adj_or_per_sd),
    direction_preserved_in_all_replicates = direction_preserved,
    n_replicates_significant = n_sig,
    relative_magnitude_change_mean = mean(rel_change),
    relative_magnitude_change_min = min(rel_change),
    relative_magnitude_change_max = max(rel_change),
    verdict = verdict,
    ## Guardrails carried on the row, not left to the writer.
    attenuation_is_not_bias = TRUE,
    interpretation = paste0(
        "Attenuation of the caudate MAGNITUDE after n-matching does not mean ",
        "the caudate estimate is biased. Module 08 tier 3 shows donor count is ",
        "a plausible major contributor to the larger magnitude. This stage asks ",
        "only whether the SCZ axis CONTRAST holds at n=", target_n, "."),
    caudate_remains_batch_confounded = TRUE,
    citable = !allow_unaccepted)

write_atomic(res, file.path(out_dir,
    paste0("scz-axis-caudate-downsampling-", opts$cohort, ".tsv")))
write_atomic(summ, file.path(out_dir,
    paste0("scz-axis-caudate-downsampling-summary-", opts$cohort, ".tsv")))

print(res[, .(comparison, cell, n_donors, n_linked, n_background,
              or_per_sd = round(adj_or_per_sd, 3),
              ci = sprintf("%.3f-%.3f", adj_or_lower, adj_or_upper),
              q = signif(adj_q, 2), significant, direction)])
cat(sprintf("\n[16] direction preserved in all %d replicates: %s\n",
            nrow(reps), direction_preserved))
cat(sprintf("[16] %d of %d replicates FDR-significant at alpha=%.2f\n",
            n_sig, nrow(reps), alpha))
cat(sprintf("[16] |log-odds| change vs full arm: mean %+.1f%% (range %+.1f%% to %+.1f%%)\n",
            100 * mean(rel_change), 100 * min(rel_change), 100 * max(rel_change)))
cat(sprintf("[16] verdict: %s\n", verdict))
cat("  Attenuation of the caudate magnitude is NOT bias in the caudate estimate.\n")
if (allow_unaccepted) {
    message("[16] --allow-unlocked: outputs are marked citable=FALSE.")
}
