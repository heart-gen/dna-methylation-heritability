#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- matched-measurability sensitivity ####
##
## Usage:
##   Rscript _h/08_matched_measurability.R --run-id rra-AA-caudate-20260906
##
## SECONDARY EVIDENCE ONLY, AND NON-GATING. Output goes to its own file and is
## never read by 03_apply_gates.R::survives(). It cannot rescue a failed primary
## and cannot break a surviving one. Config: sensitivities.matched_measurability
## (which realizes the long-declared matched_high_low_comparison key).
##
## WHY THIS EXISTS. The primary model handles technical measurability by
## REGRESSION ADJUSTMENT on gc_content, wgbs_coverage and vmr_length. That
## assumes a functional form and, worse, extrapolates across the whole covariate
## range -- including regions of it where one group has almost no support. In
## caudate that is not hypothetical: brain region is perfectly confounded with
## sequencing batch in the AANRI phase 1 delivery, and the GC-coverage
## relationship inverts between batches (see the 2026-09-06 config amendment).
## Adjustment absorbs most of that, but a proxy for a perfectly confounded
## variable is not a fix.
##
## Matching asks the cleaner question. Instead of "is the association non-zero
## after linearly adjusting for measurability", it asks: are these loci
## different from loci that were EQUALLY MEASURABLE? Two transposes are fitted,
## both declared in config:
##
##   annotation_vs_matched -- cases are loci overlapping the annotation,
##     controls are non-overlapping loci matched on GC, coverage and length.
##     The contrast is on the predictor. This is the literal form of the
##     question: are LINE/L1 loci different from equally measurable genomic loci?
##
##   score_high_vs_low -- cases are top-tertile predictor loci, controls are
##     matched bottom-tertile loci. The contrast is on the annotation fraction.
##     This is the transpose that lines up with the primary model, so a
##     disagreement between it and the primary is informative about the
##     functional form rather than about the matching.
##
## Balance is written beside the estimates. A matched analysis whose balance
## table is not reported is not interpretable, so the two files are produced in
## the same pass and neither is optional.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "04_repeat_repressive_architecture"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
feat <- fread(file.path(run_dir, "results", "vmr-features.tsv"))
annot <- load_config("repeat_annotations")

SPEC <- annot$sensitivities$matched_measurability
if (is.null(SPEC)) {
    stop("config sensitivities.matched_measurability is absent; this script is ",
         "driven by config and will not invent its own matching parameters")
}
MATCH_ON  <- unlist(SPEC$match_on)
CALIPER   <- as.numeric(SPEC$caliper_sd %||% 0.2)
RATIO     <- as.integer(SPEC$ratio %||% 1L)
SEED      <- as.integer(SPEC$seed %||% 1L)
MIN_PAIRS <- as.integer(SPEC$min_pairs %||% 100L)
CONTRASTS <- unlist(SPEC$contrasts)
if (!identical(RATIO, 1L)) {
    stop("only 1:1 matching is implemented; config asks for ratio ", RATIO)
}

PREDICTOR <- annot$primary_model$predictor %||% "local_snp_contribution_score_z"
OUTCOMES  <- c(unlist(annot$multiple_testing$family),
               names(annot$multiple_testing$outside_family))

missing <- setdiff(c(MATCH_ON, PREDICTOR, OUTCOMES), names(feat))
if (length(missing)) {
    stop("Absent from the feature table: ", paste(missing, collapse = ", "))
}

#' Standardized mean difference, the balance statistic. Reported before and
#' after matching for every match variable; a matched analysis is only as good
#' as this table.
smd <- function(a, b) {
    sd_pool <- sqrt((stats::var(a, na.rm = TRUE) + stats::var(b, na.rm = TRUE)) / 2)
    if (!is.finite(sd_pool) || sd_pool == 0) return(NA_real_)
    (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sd_pool
}

#' Greedy 1:1 nearest-neighbour matching without replacement.
#'
#' Distance is Mahalanobis on the STANDARDIZED match variables, so no variable
#' dominates through its units, and the caliper is an additional per-variable
#' constraint in the same standardized units: every matched pair must agree to
#' within `CALIPER` pooled SD on EVERY match variable, not merely on the
#' summed distance. A case with no admissible control inside the caliper is
#' dropped, and the number dropped is reported -- dropping cases silently is how
#' a matched analysis quietly becomes a different analysis.
#'
#' Greedy order is randomized under the config seed. Greedy matching is
#' order-dependent; MatchIt/optmatch (which would offer optimal matching) are
#' not installed in the project environment, so the order is fixed by seed and
#' declared rather than left to row order in the feature table.
match_pairs <- function(d, is_case) {
    z <- scale(as.matrix(d[, ..MATCH_ON]))
    ok <- stats::complete.cases(z)
    idx_case <- which(is_case & ok)
    idx_ctrl <- which(!is_case & ok)
    if (length(idx_case) == 0 || length(idx_ctrl) == 0) return(NULL)

    S <- stats::cov(z[ok, , drop = FALSE])
    Sinv <- tryCatch(solve(S), error = function(e) diag(ncol(z)))

    set.seed(SEED)
    order_case <- sample(idx_case)
    avail <- rep(TRUE, length(idx_ctrl))
    zc <- z[idx_ctrl, , drop = FALSE]

    matched_case <- integer(0)
    matched_ctrl <- integer(0)
    for (i in order_case) {
        if (!any(avail)) break
        diff <- sweep(zc, 2, z[i, ], "-")
        ## Per-variable caliper first: it is the substantive constraint, and it
        ## also shrinks the set the distance has to be computed over.
        inside <- avail & matrixStats_rowAllAbsLe(diff, CALIPER)
        if (!any(inside)) next
        cand <- which(inside)
        dm <- rowSums((diff[cand, , drop = FALSE] %*% Sinv) *
                      diff[cand, , drop = FALSE])
        pick <- cand[which.min(dm)]
        avail[pick] <- FALSE
        matched_case <- c(matched_case, i)
        matched_ctrl <- c(matched_ctrl, idx_ctrl[pick])
    }
    if (length(matched_case) == 0) return(NULL)
    list(case = matched_case, control = matched_ctrl,
         n_case_eligible = length(idx_case),
         n_ctrl_eligible = length(idx_ctrl))
}

## Small helper kept local rather than pulling in matrixStats for one call.
matrixStats_rowAllAbsLe <- function(m, k) {
    rowSums(abs(m) > k, na.rm = FALSE) == 0L
}

#' Paired inference on one matched design.
#'
#' The pairing is the point: an unpaired test on matched data throws away the
#' matching. Reported as a paired t on the differences, with a Wilcoxon signed
#' rank beside it because the contrast variables (overlap fractions especially)
#' are not close to normal. Neither is corrected -- this arm is outside the BH
#' family and is secondary evidence.
paired_result <- function(v_case, v_ctrl, contrast, outcome, contrast_var,
                          info, n_dropped) {
    delta <- v_case - v_ctrl
    delta <- delta[is.finite(delta)]
    n <- length(delta)
    if (n < MIN_PAIRS) {
        return(data.table(contrast = contrast, outcome = outcome,
                          contrast_variable = contrast_var,
                          n_pairs = n, n_cases_eligible = info$n_case_eligible,
                          n_controls_eligible = info$n_ctrl_eligible,
                          n_cases_unmatched = n_dropped,
                          estimate = NA_real_, se = NA_real_, t = NA_real_,
                          p = NA_real_, p_wilcoxon = NA_real_,
                          note = paste0("fewer than ", MIN_PAIRS, " matched pairs")))
    }
    tt <- stats::t.test(delta)
    wt <- suppressWarnings(stats::wilcox.test(delta))
    data.table(contrast = contrast, outcome = outcome,
               contrast_variable = contrast_var,
               n_pairs = n, n_cases_eligible = info$n_case_eligible,
               n_controls_eligible = info$n_ctrl_eligible,
               n_cases_unmatched = n_dropped,
               estimate = unname(tt$estimate), se = tt$stderr,
               t = unname(tt$statistic), p = tt$p.value,
               p_wilcoxon = wt$p.value, note = NA_character_)
}

balance_rows <- function(d, is_case, m, contrast, outcome) {
    rbindlist(lapply(MATCH_ON, function(v) {
        data.table(contrast = contrast, outcome = outcome, match_variable = v,
                   smd_before = smd(d[[v]][is_case], d[[v]][!is_case]),
                   smd_after  = smd(d[[v]][m$case], d[[v]][m$control]),
                   mean_case_after = mean(d[[v]][m$case], na.rm = TRUE),
                   mean_control_after = mean(d[[v]][m$control], na.rm = TRUE))
    }))
}

## ---------------------------------------------------------------- contrasts

run_annotation_vs_matched <- function(outcome) {
    any_col <- sub("_frac$", "_any", outcome)
    keep <- stats::complete.cases(feat[, c(any_col, PREDICTOR, MATCH_ON),
                                      with = FALSE])
    d <- feat[keep]
    is_case <- d[[any_col]] == 1
    if (sum(is_case) < MIN_PAIRS || sum(!is_case) < MIN_PAIRS) return(NULL)
    m <- match_pairs(d, is_case)
    if (is.null(m)) return(NULL)
    list(res = paired_result(d[[PREDICTOR]][m$case], d[[PREDICTOR]][m$control],
                             "annotation_vs_matched", outcome, PREDICTOR, m,
                             m$n_case_eligible - length(m$case)),
         bal = balance_rows(d, is_case, m, "annotation_vs_matched", outcome))
}

run_score_high_vs_low <- function(outcome) {
    keep <- stats::complete.cases(feat[, c(outcome, PREDICTOR, MATCH_ON),
                                      with = FALSE])
    d <- feat[keep]
    ## Tertiles of the predictor within this outcome's complete cases, so the
    ## split is defined on exactly the loci the contrast is fitted on.
    cuts <- stats::quantile(d[[PREDICTOR]], c(1/3, 2/3), na.rm = TRUE)
    grp <- ifelse(d[[PREDICTOR]] >= cuts[2], "high",
                  ifelse(d[[PREDICTOR]] <= cuts[1], "low", NA))
    d <- d[!is.na(grp)]
    is_case <- grp[!is.na(grp)] == "high"
    if (sum(is_case) < MIN_PAIRS || sum(!is_case) < MIN_PAIRS) return(NULL)
    m <- match_pairs(d, is_case)
    if (is.null(m)) return(NULL)
    list(res = paired_result(d[[outcome]][m$case], d[[outcome]][m$control],
                             "score_high_vs_low", outcome, outcome, m,
                             m$n_case_eligible - length(m$case)),
         bal = balance_rows(d, is_case, m, "score_high_vs_low", outcome))
}

RUNNERS <- list(annotation_vs_matched = run_annotation_vs_matched,
                score_high_vs_low = run_score_high_vs_low)
unknown <- setdiff(CONTRASTS, names(RUNNERS))
if (length(unknown)) {
    stop("config declares contrast(s) with no implementation here: ",
         paste(unknown, collapse = ", "))
}

out <- lapply(CONTRASTS, function(cn) {
    lapply(OUTCOMES, function(o) RUNNERS[[cn]](o))
})
out <- unlist(out, recursive = FALSE)
out <- out[!vapply(out, is.null, logical(1))]
if (length(out) == 0) stop("No matched contrast could be fitted for this run")

results <- rbindlist(lapply(out, `[[`, "res"))
balance <- rbindlist(lapply(out, `[[`, "bal"))

meta <- list(region = feat$region[1], population = feat$population[1],
             vmr_set_id = feat$vmr_set_id[1])
for (nm in names(meta)) {
    set(results, j = nm, value = meta[[nm]])
    set(balance, j = nm, value = meta[[nm]])
}
## Carried in the table itself, not only in this header: a reader of the TSV
## alone must be able to tell that these rows are non-gating.
results[, `:=`(gating = FALSE,
               match_on = paste(MATCH_ON, collapse = ","),
               caliper_sd = CALIPER, seed = SEED)]

write_atomic(results,
             file.path(run_dir, "results", "descriptive-matched-measurability.tsv"))
write_atomic(balance,
             file.path(run_dir, "results", "descriptive-matched-balance.tsv"))

print(results[outcome == "line_l1_frac"])
worst <- balance[which.max(abs(smd_after))]
message("[04] matched-measurability: ", nrow(results), " contrasts; ",
        "largest post-match |SMD| = ", signif(abs(worst$smd_after), 3),
        " (", worst$match_variable, ", ", worst$contrast, "/", worst$outcome, ")")
