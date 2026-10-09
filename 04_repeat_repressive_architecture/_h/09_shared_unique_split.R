#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- shared vs region-unique VMR split ####
##
## Usage:
##   Rscript _h/09_shared_unique_split.R --cohort AA [--allow-unlocked]
##
## A NON-GATING sensitivity over the three accepted cells, sealed as its own
## cross-region run (`rra-{cohort}-crossregion-{date}`). It asks one question
## about the primary fits that Figure 3 reports: is each compartment association
## a property of the region's whole VMR set, or is it carried by the VMRs that
## only that region called?
##
## Why it exists. The 2026-09-23 audit refit the primary model on the split that
## _h/07_qc_covariate_attribution.R declares (pass `set_overlap`) and found
## LINE/L1 and H3K9me3 null on VMRs shared by all three regions while quiescent
## and accessible held at about half strength. That refit was scratch: unsealed,
## run on hand-picked tables, not citable. The Figure 3 sentence cannot rest on
## it, and it cannot be written without it (TASKS.md A5). This stage makes it a
## run.
##
## What it is NOT. It changes no gate, no claim and no q-value: the accepted
## cells, their interpretation-claims.tsv and their decision tokens are read
## only. A subset of a region's VMRs is not a locked analysis set, so nothing
## here enters 03_apply_gates.R's conjunction. Every row carries
## gating = FALSE.
##
## The split, unchanged from _h/07 so the two can be compared:
##   shared_with_both  the VMR overlaps (>= 1 bp) a VMR in EACH other region
##   region_unique     every other VMR of the region
## Overlap is physical, on the accepted catalogs. It is a property of which loci
## were CALLED in each region, and region is confounded with sequencing batch for
## caudate (AGENTS.md 8.1), so "region-unique" in caudate also means
## "batch-3-unique". Caudate LINE/L1 stays set aside exactly as in the primary.
##
## The model is the primary model of _h/02_test_association.R, rebuilt from the
## same config and checked, before any subset is fitted, to reproduce every
## sealed primary estimate to 1e-6. If it does not, the stage stops: a subset
## result from a different model than the one Figure 3 reports answers a
## different question.
##
## Inference. Each subset estimate carries its HC3 SE (the primary's) and a
## delete-one-chromosome weighted block-jackknife SE (Busing et al. 1999; the
## construction 00_shared/axis_inference.R uses). The unique-minus-shared
## difference is tested with the JOINT jackknife -- chromosome k removed from
## both subsets at once -- because the two subsets are disjoint in loci but not
## independent: nearby VMRs share sequence context and every VMR is estimated in
## the same donors. The HC3 independence SE is reported beside it for reference
## and is not the test.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(sandwich)
    library(lmtest)
    library(GenomicRanges)
})

MODULE <- "04_repeat_repressive_architecture"
MODULE_TAG <- "rra"
PRED <- "local_snp_contribution_score_z"
GROUPS <- c("shared_with_both", "region_unique")
MIN_N <- 200L        # _h/07's floor for a subset fit
REPRO_TOL <- 1e-6

opts <- parse_v2_args(require = "cohort")
allow_unlocked <- isTRUE(opts$allow_unlocked)
cohort <- opts$cohort
if (!identical(parse_cell(cohort)$cell_kind, "arm")) {
    stop("The split runs on a discovery arm's catalogs, not the cell '", cohort, "'")
}

annot <- load_config("repeat_annotations")
assert_locked(list(repeat_annotations = annot), allow_unlocked = allow_unlocked)
regions <- c("caudate", "dlpfc", "hippocampus")

## ------------------------------------------------------------ the three cells
up <- lapply(regions, function(re)
    require_accepted_upstream(MODULE, cohort, re, allow_unaccepted = allow_unlocked))
names(up) <- regions
if (any(is.na(vapply(up, function(u) as.character(u$run_id), character(1))))) {
    stop("A smoke split needs accepted cells to read; there is no smoke input path.")
}
cell_dir <- function(re) file.path(repo_root(), MODULE, "_m", "runs", up[[re]]$run_id)

## The config the cells were fitted under must be the config this stage reads,
## or the reproduction guard below is the only thing standing between a changed
## adjustment set and a silently different model. Refuse early and say why.
cfg_sha <- attr(annot, "config_sha256")
for (re in regions) {
    m <- fread(file.path(cell_dir(re), "manifest.tsv"), colClasses = "character")
    rec <- m$value[m$field == "config_repeat_annotations_sha256"]
    if (!length(rec) || !identical(rec[1], cfg_sha)) {
        stop(up[[re]]$run_id, " was fitted under config_repeat_annotations_sha256 ",
             rec[1], "; config/repeat_annotations.yml is now ", cfg_sha,
             ". Refit the cells, or run this stage from the commit they record.")
    }
}

feats <- lapply(regions, function(re)
    fread(file.path(cell_dir(re), "results", "vmr-features.tsv")))
names(feats) <- regions
sealed <- lapply(regions, function(re)
    fread(file.path(cell_dir(re), "results", "association-results.tsv"))[
        analysis_set == "primary" & predictor == PRED])
names(sealed) <- regions

## ---------------------------------------------------------- the primary model
## Mirrors _h/02_test_association.R: config-declared covariates mapped to their
## realized terms, minus any the cell recorded as dropped for being constant.
REALIZED_TERM <- list(
    vmr_length = "log(vmr_length)", cpg_count = "cpg_count",
    cpg_density = "cpg_density", gc_content = "gc_content",
    mean_methylation = "mean_methylation",
    methylation_variance = "methylation_variance",
    wgbs_coverage = "wgbs_coverage", tested_snp_count = "tested_snp_count",
    snp_proximity = "snp_proximal_frac", mappability = "mappability",
    segdup_overlap = "segdup_frac",
    problematic_region_overlap = "problematic_frac",
    broad_genomic_annotation = "factor(broad_genomic_annotation)",
    cell_composition_pcs = "cell_composition_r2"
)
declared <- unlist(annot$covariates)
unknown <- setdiff(declared, names(REALIZED_TERM))
if (length(unknown)) stop("Covariate(s) with no model term: ", paste(unknown, collapse = ", "))
term_column <- function(t) gsub("^(log|factor)\\(|\\)$", "", t)

covariates_for <- function(re) {
    terms <- unlist(REALIZED_TERM[declared])
    dropped_f <- file.path(cell_dir(re), "results", "dropped-covariates.tsv")
    if (file.exists(dropped_f)) {
        dropped <- fread(dropped_f)$covariate
        terms <- terms[!names(terms) %in% dropped]
    }
    unname(terms)
}

BH_FAMILY <- unlist(annot$multiple_testing$family)
CONTROLS <- names(annot$multiple_testing$outside_family)
OUTCOMES <- c(BH_FAMILY, CONTROLS)          # the fraction scale Figure 3 draws
role_of <- function(o) {
    if (o %in% BH_FAMILY) return("bh_family")
    annot$multiple_testing$outside_family[[o]]$role %||% "control"
}
set_aside <- annot$interpretation$technically_confounded_regions

## A subset can make a covariate constant (a factor level that never occurs, a
## segdup fraction that is zero throughout). Such a term is dropped FOR THAT FIT
## and recorded, rather than left to alias and break the sandwich.
usable_terms <- function(d, terms) {
    keep <- vapply(terms, function(t) {
        v <- d[[term_column(t)]]
        length(unique(v[!is.na(v)])) >= 2L
    }, logical(1))
    list(keep = terms[keep], dropped = terms[!keep])
}

fit_glm <- function(d, outcome, terms) {
    form <- stats::as.formula(paste(outcome, "~", PRED, "+", paste(terms, collapse = " + ")))
    fam <- stats::quasibinomial(link = annot$primary_model$link %||% "logit")
    tryCatch(suppressWarnings(stats::glm(form, data = d, family = fam)),
             error = function(e) NULL)
}

coef_pred <- function(fit) {
    if (is.null(fit)) return(NA_real_)
    b <- stats::coef(fit)[PRED]
    if (is.na(b)) NA_real_ else unname(b)
}

fit_hc3 <- function(d, outcome, terms) {
    fit <- fit_glm(d, outcome, terms)
    if (is.null(fit) || is.na(stats::coef(fit)[PRED])) {
        return(list(estimate = NA_real_, se = NA_real_, p = NA_real_, n_fitted = 0L))
    }
    ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(
        fit, type = annot$primary_model$vcov %||% "HC3"))
    list(estimate = ct[PRED, 1], se = ct[PRED, 2], p = ct[PRED, 4],
         n_fitted = as.integer(stats::nobs(fit)))
}

## Weighted delete-one-chromosome jackknife over the rows of `d`. `theta_fn`
## maps a row subset to the statistic, so the same code serves one subset's
## coefficient and the joint unique-minus-shared difference.
block_jackknife <- function(chrom, theta_fn, full) {
    blocks <- sort(unique(chrom))
    if (length(blocks) < 2L || !is.finite(full)) return(NA_real_)
    hj <- length(chrom) / vapply(blocks, function(b) sum(chrom == b), numeric(1))
    theta <- vapply(blocks, function(b) theta_fn(chrom != b), numeric(1))
    ok <- is.finite(theta)
    if (sum(ok) < 2L) return(NA_real_)
    pseudo <- hj[ok] * full - (hj[ok] - 1) * theta[ok]
    sqrt(sum((pseudo - mean(pseudo))^2 / (hj[ok] - 1)) / sum(ok))
}

## ------------------------------------------------------- reproduction guard
repro <- rbindlist(lapply(regions, function(re) {
    d <- feats[[re]]
    terms <- covariates_for(re)
    rbindlist(lapply(OUTCOMES, function(o) {
        r <- fit_hc3(d, o, terms)
        s <- sealed[[re]][outcome == o]
        if (nrow(s) != 1L) stop("No unique sealed primary row for ", re, " x ", o)
        data.table(region = re, outcome = o,
                   estimate_refit = r$estimate, estimate_sealed = s$estimate,
                   se_refit = r$se, se_sealed = s$se,
                   n_fitted_refit = r$n_fitted, n_fitted_sealed = s$n_fitted)
    }))
}))
repro[, `:=`(abs_diff_estimate = abs(estimate_refit - estimate_sealed),
             abs_diff_se = abs(se_refit - se_sealed))]
bad <- repro[!(abs_diff_estimate <= REPRO_TOL & abs_diff_se <= REPRO_TOL &
               n_fitted_refit == n_fitted_sealed)]
if (nrow(bad)) {
    print(bad)
    stop(nrow(bad), " primary fit(s) do not reproduce the sealed cells. The split ",
         "would be fitted on a different model from the one Figure 3 reports.")
}
message("[04] reproduction guard: ", nrow(repro), " primary fits reproduce to <= ",
        REPRO_TOL)

## ----------------------------------------------------------------- the split
as_gr <- function(d) GRanges(d$chrom, IRanges(d$start, d$end))
membership <- rbindlist(lapply(regions, function(re) {
    d <- feats[[re]]
    others <- setdiff(regions, re)
    hit <- vapply(others, function(o) overlapsAny(as_gr(d), as_gr(feats[[o]])),
                  logical(nrow(d)))
    data.table(region = re, vmr_id = d$vmr_id, chrom = d$chrom,
               overlaps_other_regions = rowSums(hit),
               vmr_group = fifelse(rowSums(hit) == length(others),
                                   "shared_with_both", "region_unique"))
}))

## What differs between the two groups, so a reader can see whether a gap in the
## slope travels with a gap in the covariates the model adjusts for.
DESCRIBE <- c(PRED, "gc_content", "mappability", "wgbs_coverage", "cpg_density",
              "mean_methylation", "segdup_frac", OUTCOMES)
composition <- rbindlist(lapply(regions, function(re) {
    d <- merge(feats[[re]], membership[region == re, .(vmr_id, vmr_group)], by = "vmr_id")
    d[, vmr_length := end - start + 1L]
    rbindlist(lapply(GROUPS, function(g) {
        s <- d[vmr_group == g]
        cols <- intersect(c(DESCRIBE, "vmr_length"), names(s))
        data.table(region = re, vmr_group = g, n_vmrs = nrow(s),
                   variable = cols,
                   mean = vapply(cols, function(k) mean(s[[k]], na.rm = TRUE), numeric(1)),
                   sd = vapply(cols, function(k) stats::sd(s[[k]], na.rm = TRUE), numeric(1)))
    }))
}))

## ------------------------------------------------------------- the subset fits
assoc <- list(); diffs <- list()
for (re in regions) {
    d <- merge(feats[[re]], membership[region == re, .(vmr_id, vmr_group)], by = "vmr_id")
    terms_all <- covariates_for(re)
    for (o in OUTCOMES) {
        sub_est <- list()
        for (g in c("all", GROUPS)) {
            s <- if (g == "all") d else d[vmr_group == g]
            if (nrow(s) < MIN_N) {
                assoc[[length(assoc) + 1]] <- data.table(
                    region = re, outcome = o, vmr_group = g, n = nrow(s),
                    note = paste("fewer than", MIN_N, "loci"))
                next
            }
            ut <- usable_terms(s, terms_all)
            r <- fit_hc3(s, o, ut$keep)
            ## Chromatin outcomes are NA off the hg19 liftover; jackknife over
            ## the rows the model can use, so a block's weight is its real n.
            s_ok <- s[is.finite(get(o))]
            se_jk <- block_jackknife(
                s_ok$chrom,
                function(keep) coef_pred(fit_glm(s_ok[keep], o, ut$keep)),
                r$estimate)
            sub_est[[g]] <- list(estimate = r$estimate, se = r$se, rows = s_ok)
            assoc[[length(assoc) + 1]] <- data.table(
                region = re, outcome = o, vmr_group = g, n = nrow(s),
                n_fitted = r$n_fitted, estimate = r$estimate,
                se_hc3 = r$se, p_hc3 = r$p, se_jackknife = se_jk,
                p_jackknife = 2 * stats::pnorm(-abs(r$estimate / se_jk)),
                terms_dropped_in_subset = paste(ut$dropped, collapse = ";"),
                note = NA_character_)
        }
        if (!all(GROUPS %in% names(sub_est))) next
        u <- sub_est$region_unique; s <- sub_est$shared_with_both
        delta <- u$estimate - s$estimate
        both <- rbind(u$rows, s$rows)
        ut_u <- usable_terms(u$rows, terms_all)$keep
        ut_s <- usable_terms(s$rows, terms_all)$keep
        se_joint <- block_jackknife(
            both$chrom,
            function(keep) {
                kb <- both[keep]
                coef_pred(fit_glm(kb[vmr_group == "region_unique"], o, ut_u)) -
                    coef_pred(fit_glm(kb[vmr_group == "shared_with_both"], o, ut_s))
            },
            delta)
        diffs[[length(diffs) + 1]] <- data.table(
            region = re, outcome = o,
            estimate_region_unique = u$estimate, estimate_shared = s$estimate,
            difference_unique_minus_shared = delta,
            se_joint_jackknife = se_joint,
            z = delta / se_joint,
            p = 2 * stats::pnorm(-abs(delta / se_joint)),
            se_if_independent_hc3 = sqrt(u$se^2 + s$se^2))
    }
    message("[04] ", re, " split fitted")
}
assoc <- rbindlist(assoc, fill = TRUE)
diffs <- rbindlist(diffs, fill = TRUE)

## Labels every row carries, so the TSV read alone cannot be mistaken for a gate.
label <- function(dt) {
    dt[, outcome_role := vapply(outcome, role_of, character(1))]
    dt[, set_aside_from_claim := mapply(function(o, re) re %in% unlist(set_aside[[o]]),
                                        outcome, region)]
    dt[, `:=`(predictor = PRED, gating = FALSE,
              sensitivity = "shared_vs_region_unique_vmrs",
              population = cohort,
              source_run_id = vapply(region, function(re) up[[re]]$run_id, character(1)),
              vmr_set_id = vapply(region, function(re) up[[re]]$vmr_set_id, character(1)))]
    dt[]
}
assoc <- label(assoc)
diffs <- label(diffs)

## ------------------------------------------------------------------ the run
if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked)")
}
run <- new_run(
    module = MODULE_TAG, cohort = cohort, region = "crossregion",
    module_root = file.path(repo_root(), MODULE), run_id = opts$run_id,
    upstream = stats::setNames(lapply(regions, function(re) up[[re]]$run_id),
                               paste0("repeat_architecture_", regions)),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_repeat_annotations_sha256 = cfg_sha,
        sensitivity = "shared_vs_region_unique_vmrs",
        gating = "FALSE",
        split_definition = "shared_with_both = overlaps >=1 bp a VMR in each other region",
        min_subset_n = MIN_N,
        reproduction_tolerance = REPRO_TOL
    )
)
res_dir <- file.path(run$dir, "results")
dir.create(res_dir, showWarnings = FALSE)
code_dir <- file.path(run$dir, "code", "_h")
dir.create(code_dir, recursive = TRUE)
invisible(file.copy(file.path(repo_root(), MODULE, "_h", "09_shared_unique_split.R"), code_dir))

write_atomic(membership, file.path(res_dir, "split-membership.tsv"))
write_atomic(composition, file.path(res_dir, "split-composition.tsv"))
write_atomic(repro, file.path(res_dir, "split-reproduction-check.tsv"))
write_atomic(assoc, file.path(res_dir, "split-association.tsv"))
write_atomic(diffs, file.path(res_dir, "split-difference.tsv"))
writeLines(c(
    "Interpretation constraints carried by this run:",
    "  - NON-GATING. No gate, claim, q-value or decision token of the accepted",
    "    cells changes. The split is a description of where each primary",
    "    association lives, not a new test of it.",
    "  - 'Shared' and 'region-unique' describe which loci were CALLED. Region is",
    "    confounded with sequencing batch for caudate (AGENTS.md 8.1), so a",
    "    caudate-unique VMR is also batch-3-unique.",
    "  - Two readings of a gap fit these data and this run cannot separate them:",
    "    the association is genuinely region-restricted, or region-unique loci are",
    "    enriched for loci whose calling, coverage or score is sensitive to local",
    "    sequence. Neither licenses a cell-type or retrotransposition claim.",
    "  - The difference test is the JOINT chromosome jackknife; the HC3",
    "    independence SE is a reference and not the test.",
    "  - Caudate LINE/L1 remains set aside from the claim (set_aside_from_claim)."
), file.path(res_dir, "interpretation-constraints.txt"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "session-info.txt"))

if (!allow_unlocked) close_run(run)
message("[04] split run ", run$run_id, if (allow_unlocked) " (smoke, not sealed)" else " sealed")
print(diffs[outcome %in% BH_FAMILY | outcome %in% c("accessible_frac", "h3k27me3_frac"),
            .(region, outcome, estimate_shared, estimate_region_unique,
              difference_unique_minus_shared, se_joint_jackknife, p)], digits = 3)
cat(run$run_id, "\n", sep = "")
