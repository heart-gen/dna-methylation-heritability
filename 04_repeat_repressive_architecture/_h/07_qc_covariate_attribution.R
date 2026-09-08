#!/usr/bin/env Rscript
##
## Post-hoc QC: WHICH covariate carries the attenuation, and is it the same
## covariate in every region?
##
##   Rscript _h/07_qc_covariate_attribution.R
##   Rscript _h/07_qc_covariate_attribution.R --outcome line_l1_frac
##
## _h/06_qc_adjustment_ladder.R decomposes the v1 -> v2 attenuation into
## cumulative BLOCKS. That is the right resolution for the confounder vs
## mediator judgement, but it cannot say which member of a block does the work,
## and it cannot show that a covariate acts in OPPOSITE directions in different
## regions. Both matter for LINE/L1: the caudate estimate falls 0.746 -> 0.038
## between v1_like and the prespecified primary while DLPFC and hippocampus hold
## at ~0.3, and the block ladder attributes that to plus_cpg + plus_seq without
## naming the variable.
##
## Four passes, each answering one question:
##   attribution  -- add each covariate singly to v1_like, and drop each singly
##                   from the prespecified primary. A covariate that attenuates
##                   in one region and SUPPRESSES in another is visible here and
##                   nowhere in the block ladder.
##   correlation  -- Spearman correlations of the predictor and the outcome with
##                   every covariate, by region. Says whether a region-specific
##                   attenuation is a property of the PREDICTOR rather than of
##                   the outcome annotation.
##   coverage     -- refit over a sweep of minimum wgbs_coverage thresholds.
##                   Tests directly whether a region-specific entanglement is a
##                   marginal-coverage artifact that a stricter Module 01 filter
##                   would remove. It is an approximation: the real filter is
##                   per-CpG in 01_vmr_catalog/_h/00_prepare.R and re-running it
##                   would recall the VMRs, whereas this only drops whole VMRs
##                   whose mean depth is low. It is therefore a screen, and a
##                   negative result here is the informative one.
##   set_overlap  -- split each region's VMRs into those shared with both other
##                   regions and those unique to it. Says whether an anomaly is
##                   carried by the region-specific part of the VMR set or is a
##                   whole-region property.
##
## This is post-hoc analysis OF sealed runs, not part of them. It writes to
## _m/qc/{run_id}/, which is gitignored and regenerable.
##
## Read the attribution table as a decomposition, NOT as a menu. Choosing a
## per-region adjustment set from these numbers is exactly the move the
## 2026-09-02 config amendment warns against; the covariates are locked in
## config/repeat_annotations.yml and this script does not fit anything that is
## allowed to become a primary estimate.

suppressPackageStartupMessages({
    library(data.table); library(sandwich); library(lmtest)
    library(GenomicRanges)
})
source(file.path(Sys.getenv("V2_REPO", "."), "00_shared", "load.R"))

MODULE <- "04_repeat_repressive_architecture"
PRED   <- "local_snp_contribution_score_z"

## Must match _h/02_test_association.R after the 2026-09-02 amendment, including
## log(vmr_length). These are the nine terms of the prespecified primary.
BASE <- c("log(vmr_length)", "tested_snp_count")
PRESPECIFIED <- c(BASE, "cpg_count", "cpg_density", "gc_content",
                  "snp_proximal_frac", "mappability", "segdup_frac",
                  "problematic_frac")
## Descriptive covariates: not in the primary formula, but the correlation and
## coverage passes need them to explain WHY a primary covariate behaves oddly.
DESCRIPTIVE <- c("mean_methylation", "methylation_variance", "wgbs_coverage",
                 "cell_composition_r2")

COV_THRESHOLDS <- c(0, 8, 10, 12, 15, 20)
MIN_N <- 200L

fit_terms <- function(d, outcome, terms) {
    form <- as.formula(paste(outcome, "~", PRED, "+",
                             paste(terms, collapse = " + ")))
    fit <- tryCatch(glm(form, data = d, family = quasibinomial(link = "logit")),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ct <- tryCatch(coeftest(fit, vcov. = vcovHC(fit, type = "HC3")),
                   error = function(e) NULL)
    if (is.null(ct) || !PRED %in% rownames(ct)) return(NULL)
    data.table(n_fitted = nobs(fit), estimate = ct[PRED, 1],
               se = ct[PRED, 2], p = ct[PRED, 4])
}

tag <- function(dt, ...) if (is.null(dt)) NULL else cbind(data.table(...), dt)

## ------------------------------------------------------------- attribution
## Two directions, because neither alone identifies the operative covariate: a
## covariate can look harmless added to a thin model yet be the only thing
## holding the full model down, and vice versa.
attribution_one <- function(d, run_id, region, outcome) {
    singles <- setdiff(PRESPECIFIED, BASE)
    rbind(
        tag(fit_terms(d, outcome, BASE), run_id = run_id, region = region,
            outcome = outcome, pass = "attribution", model = "v1_like",
            term = NA_character_),
        tag(fit_terms(d, outcome, PRESPECIFIED), run_id = run_id,
            region = region, outcome = outcome, pass = "attribution",
            model = "prespecified", term = NA_character_),
        rbindlist(lapply(singles, function(tm)
            tag(fit_terms(d, outcome, c(BASE, tm)), run_id = run_id,
                region = region, outcome = outcome, pass = "attribution",
                model = "v1_like_plus_one", term = tm))),
        rbindlist(lapply(singles, function(tm)
            tag(fit_terms(d, outcome, setdiff(PRESPECIFIED, tm)),
                run_id = run_id, region = region, outcome = outcome,
                pass = "attribution", model = "prespecified_minus_one",
                term = tm))),
        fill = TRUE)
}

## ------------------------------------------------------------- correlation
## Reported for the predictor AND the outcome. If a covariate is entangled with
## the predictor in one region only, the attenuation there is about how the
## score was constructed, not about the compartment being tested.
correlation_one <- function(d, run_id, region, outcome) {
    vars <- intersect(c(sub("^log\\(|\\)$", "", PRESPECIFIED), DESCRIPTIVE),
                      names(d))
    rbindlist(lapply(vars, function(v) data.table(
        run_id = run_id, region = region, outcome = outcome,
        pass = "correlation", model = "spearman", term = v,
        cor_with_predictor = cor(d[[v]], d[[PRED]], method = "spearman",
                                 use = "pairwise.complete.obs"),
        cor_with_outcome = cor(d[[v]], d[[outcome]], method = "spearman",
                               use = "pairwise.complete.obs"),
        cor_with_coverage = if ("wgbs_coverage" %in% names(d))
            cor(d[[v]], d$wgbs_coverage, method = "spearman",
                use = "pairwise.complete.obs") else NA_real_)))
}

## ------------------------------------------------------------- coverage
coverage_one <- function(d, run_id, region, outcome) {
    if (!"wgbs_coverage" %in% names(d)) return(NULL)
    rbindlist(lapply(COV_THRESHOLDS, function(thr) {
        s <- d[wgbs_coverage >= thr]
        if (nrow(s) < MIN_N) return(NULL)
        rbind(
            tag(fit_terms(s, outcome, BASE), run_id = run_id, region = region,
                outcome = outcome, pass = "coverage", model = "v1_like",
                term = paste0("wgbs_coverage>=", thr)),
            tag(fit_terms(s, outcome, PRESPECIFIED), run_id = run_id,
                region = region, outcome = outcome, pass = "coverage",
                model = "prespecified",
                term = paste0("wgbs_coverage>=", thr)),
            fill = TRUE)[, `:=`(
                min_coverage = thr, n_set = nrow(s),
                cor_pred_gc = cor(s[[PRED]], s$gc_content, method = "spearman",
                                  use = "pairwise.complete.obs"),
                cor_gc_coverage = cor(s$gc_content, s$wgbs_coverage,
                                      method = "spearman",
                                      use = "pairwise.complete.obs"))]
    }), fill = TRUE)
}

## ------------------------------------------------------------- set_overlap
as_gr <- function(d) GRanges(d$chrom, IRanges(d$start, d$end))

overlap_one <- function(d, others, run_id, region, outcome) {
    shared <- Reduce(`&`, lapply(others, function(o)
        overlapsAny(as_gr(d), as_gr(o))))
    d <- copy(d)[, vmr_group := fifelse(shared, "shared_with_both",
                                        "region_unique")]
    rbindlist(lapply(c("shared_with_both", "region_unique"), function(g) {
        s <- d[vmr_group == g]
        if (nrow(s) < MIN_N) return(NULL)
        rbind(
            tag(fit_terms(s, outcome, BASE), run_id = run_id, region = region,
                outcome = outcome, pass = "set_overlap", model = "v1_like",
                term = g),
            tag(fit_terms(s, outcome, PRESPECIFIED), run_id = run_id,
                region = region, outcome = outcome, pass = "set_overlap",
                model = "prespecified", term = g),
            fill = TRUE)[, `:=`(
                n_set = nrow(s), mean_gc = mean(s$gc_content, na.rm = TRUE),
                mean_coverage = if ("wgbs_coverage" %in% names(s))
                    mean(s$wgbs_coverage, na.rm = TRUE) else NA_real_,
                cor_pred_gc = cor(s[[PRED]], s$gc_content, method = "spearman",
                                  use = "pairwise.complete.obs"),
                cor_gc_coverage = if ("wgbs_coverage" %in% names(s))
                    cor(s$gc_content, s$wgbs_coverage, method = "spearman",
                        use = "pairwise.complete.obs") else NA_real_)]
    }), fill = TRUE)
}

## ------------------------------------------------------------------- driver
args <- commandArgs(trailingOnly = TRUE)
OUTCOME <- if (length(args) >= 2 && args[1] == "--outcome") args[2] else
    "line_l1_frac"

runs_dir <- file.path(repo_root(), MODULE, "_m", "runs")
run_ids <- grep("^rra-AA-[a-z]+-[0-9]{8}$",
                list.dirs(runs_dir, recursive = FALSE, full.names = FALSE),
                value = TRUE)
if (length(run_ids) == 0) stop("no production run ids resolved under ", runs_dir)
## The set_overlap pass is cross-region by construction, so a single-region
## invocation cannot support it; the other three passes are per-region.
if (length(run_ids) < 3) message("[07] fewer than three regions; ",
                                 "set_overlap pass will be skipped")

feats <- lapply(run_ids, function(r)
    fread(file.path(runs_dir, r, "results", "vmr-features.tsv")))
names(feats) <- run_ids
regions <- vapply(run_ids, function(r)
    sub("^rra-[A-Za-z]+-([a-z]+)-.*$", "\\1", r), character(1))

res <- rbindlist(lapply(seq_along(run_ids), function(i) {
    d <- feats[[i]]; rid <- run_ids[i]; reg <- regions[i]
    if (!OUTCOME %in% names(d))
        stop("outcome ", OUTCOME, " absent from ", rid)
    out <- list(attribution_one(d, rid, reg, OUTCOME),
                correlation_one(d, rid, reg, OUTCOME),
                coverage_one(d, rid, reg, OUTCOME))
    if (length(run_ids) >= 3)
        out[[length(out) + 1]] <- overlap_one(d, feats[-i], rid, reg, OUTCOME)
    rbindlist(out, fill = TRUE)
}), fill = TRUE)

setorder(res, pass, outcome, model, term, region)

qc_dir <- file.path(repo_root(), MODULE, "_m", "qc", run_ids[1])
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(qc_dir, paste0("covariate-attribution-", OUTCOME, ".tsv"))
fwrite(res, outfile, sep = "\t")

show <- function(p, ttl, vv = c("estimate", "p")) {
    x <- res[pass == p]
    if (!nrow(x)) return(invisible(NULL))
    cat("\n===", ttl, "--", OUTCOME, "\n")
    print(dcast(x, model + term ~ region, value.var = vv), digits = 3)
}
show("attribution", "per-covariate attribution")
cat("\n=== predictor/outcome correlations --", OUTCOME, "\n")
print(dcast(res[pass == "correlation"], term ~ region,
            value.var = c("cor_with_predictor", "cor_with_outcome")), digits = 3)
show("coverage", "coverage-threshold sweep",
     c("n_set", "estimate", "p", "cor_pred_gc", "cor_gc_coverage"))
show("set_overlap", "shared vs region-unique VMRs",
     c("n_set", "estimate", "p", "cor_pred_gc", "cor_gc_coverage"))
cat("\nwrote", outfile, "\n")
