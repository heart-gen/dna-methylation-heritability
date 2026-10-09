#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- LINE/L1 subfamily resolution ####
##
## Usage:
##   Rscript _h/10_l1_subfamily.R --cohort AA [--allow-unlocked]
##
## A NON-GATING secondary breakdown of line_l1_frac over the three accepted
## cells, sealed as its own cross-region run (`rra-{cohort}-crossregion-{date}`).
## It realizes the AGENTS.md 7.4 high-value extension: separate LINE/L1 by
## subfamily age (young L1HS/L1PA vs old L1M) and by completeness (full-length vs
## fragment). PI direction 2026-10-07; definitions in config/l1_subfamilies.yml.
##
## What it is NOT. It changes no gate, claim or q-value of the accepted cells:
## their interpretation-claims.tsv and decision tokens are read only, and every
## row written here carries gating = FALSE. Overlap with a young or full-length
## element is overlap. It is not activity, expression, or retrotransposition, and
## full-length is not retrotransposition competence (no ORF check).
##
## Two guards run before anything is fitted, and either one stops the stage:
##   1. PARTITION. The subfamily classes are recomputed per VMR from the
##      subfamily asset, and their union must reproduce the sealed
##      line_l1_frac of every VMR to 1e-12. That proves the asset and the
##      overlap arithmetic are the ones the accepted cells used.
##   2. MODEL. The primary line_l1_frac fit is rebuilt from the same config and
##      must reproduce the sealed estimate and SE to 1e-6.
##
## Inference. Each estimate carries its HC3 SE (the primary's) and a
## delete-one-chromosome weighted block-jackknife SE. Each contrast is a
## difference of two slopes on the SAME VMRs, so the two are not independent;
## it is tested with the JOINT jackknife (chromosome k removed from both fits at
## once), the construction _h/09_shared_unique_split.R uses.

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
REPRO_TOL <- 1e-6
PARTITION_TOL <- 1e-12

opts <- parse_v2_args(require = "cohort")
allow_unlocked <- isTRUE(opts$allow_unlocked)
cohort <- opts$cohort
if (!identical(parse_cell(cohort)$cell_kind, "arm")) {
    stop("The subfamily stage runs on a discovery arm's cells, not the cell '", cohort, "'")
}

annot <- load_config("repeat_annotations")
l1cfg <- load_config("l1_subfamilies")
assert_locked(list(repeat_annotations = annot, l1_subfamilies = l1cfg),
              allow_unlocked = allow_unlocked)
regions <- c("caudate", "dlpfc", "hippocampus")

## ------------------------------------------------------------ the three cells
up <- lapply(regions, function(re)
    require_accepted_upstream(MODULE, cohort, re, allow_unaccepted = allow_unlocked))
names(up) <- regions
if (any(is.na(vapply(up, function(u) as.character(u$run_id), character(1))))) {
    stop("The subfamily stage reads accepted cells; there is no smoke input path.")
}
cell_dir <- function(re) file.path(repo_root(), MODULE, "_m", "runs", up[[re]]$run_id)

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
        analysis_set == "primary" & predictor == PRED & outcome == "line_l1_frac"])
names(sealed) <- regions

## ------------------------------------------------------- the subfamily asset
asset_f <- file.path(repo_root(), l1cfg$source$subfamily_asset)
if (!file.exists(asset_f)) {
    stop("Missing ", asset_f, ". Build it with ",
         "inputs/supportfiles/_h/02_build_l1_subfamily_asset.py")
}
asset_sha <- file_sha256(asset_f)
l1 <- fread(cmd = paste("zcat", shQuote(asset_f)))
l1[, full_length := as.logical(full_length)]
l1[, retains_5prime := as.logical(retains_5prime)]
## BED is 0-based half-open; annotation_io.R reads the L1 BED through
## rtracklayer::import(), which adds 1 to the start. Do the same.
l1_gr <- GRanges(l1$chrom, IRanges(l1$start + 1L, l1$end))

class_rows <- function(spec) {
    keep <- rep(TRUE, nrow(l1))
    if (!is.null(spec$age)) keep <- keep & l1$age_class == spec$age
    if (!is.null(spec$full_length)) keep <- keep & l1$full_length == isTRUE(spec$full_length)
    if (!is.null(spec$retains_5prime)) keep <- keep & l1$retains_5prime == isTRUE(spec$retains_5prime)
    keep
}
OUTCOMES <- names(l1cfg$outcomes)
CONTRASTS <- l1cfg$contrasts
class_gr <- lapply(l1cfg$outcomes, function(spec) l1_gr[class_rows(spec)])
asset_counts <- data.table(
    outcome = OUTCOMES,
    n_elements = vapply(class_gr, length, integer(1)),
    bp_reduced = vapply(class_gr, function(g)
        sum(as.numeric(width(GenomicRanges::reduce(g, ignore.strand = TRUE)))), numeric(1)))

## Summed overlap fraction, as annotation_io.R::overlap_features() computes it.
overlap_frac <- function(gr, ann_gr) {
    ann <- GenomicRanges::reduce(ann_gr, ignore.strand = TRUE)
    hits <- findOverlaps(gr, ann, ignore.strand = TRUE)
    frac <- numeric(length(gr))
    if (length(hits) > 0) {
        inter <- pintersect(gr[queryHits(hits)], ann[subjectHits(hits)],
                            ignore.strand = TRUE)
        cov <- tapply(width(inter), queryHits(hits), sum)
        frac[as.integer(names(cov))] <- as.numeric(cov)
    }
    pmin(frac / width(gr), 1)
}

## ------------------------------------------------------- guard 1: partition
partition <- rbindlist(lapply(regions, function(re) {
    d <- feats[[re]]
    gr <- GRanges(d$chrom, IRanges(d$start, d$end))
    union_frac <- overlap_frac(gr, l1_gr)
    for (o in OUTCOMES) set(d, j = o, value = overlap_frac(gr, class_gr[[o]]))
    set(d, j = "l1_union_recomputed", value = union_frac)
    feats[[re]] <<- d
    age_sum <- d$l1_young_frac + d$l1_old_frac + d$l1_intermediate_frac
    comp_sum <- d$l1_full_length_frac + d$l1_fragment_frac
    data.table(region = re, n_vmrs = nrow(d),
               max_abs_union_vs_sealed = max(abs(union_frac - d$line_l1_frac)),
               max_abs_age_sum_vs_sealed = max(abs(age_sum - d$line_l1_frac)),
               max_abs_completeness_sum_vs_sealed = max(abs(comp_sum - d$line_l1_frac)))
}))
print(partition)
if (any(partition$max_abs_union_vs_sealed > PARTITION_TOL)) {
    stop("The subfamily asset does not reproduce the sealed line_l1_frac; it is ",
         "not the annotation the accepted cells used.")
}
## The classes are disjoint SETS of elements, but two elements of different
## classes can abut or overlap by a few bp in the genome, so a class sum may
## exceed the union slightly. Record it; it is not a failure.

## ---------------------------------------------------------- the primary model
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

## Analysis sets as 02_test_association.R realizes them.
SETS <- list(
    primary = list(subset = function(d) d, extra = character(0)),
    high_mappability = list(
        subset = function(d)
            d[mappability >= annot$sensitivities$high_mappability$min_mappability],
        extra = character(0)),
    exclude_segdups = list(subset = function(d) d[segdup_frac == 0], extra = character(0)),
    adjust_cell_composition = list(subset = function(d) d, extra = "cell_composition_r2")
)
unknown_sets <- setdiff(unlist(l1cfg$analysis_sets), names(SETS))
if (length(unknown_sets)) stop("Unrealized analysis set(s): ", paste(unknown_sets, collapse = ", "))
SETS <- SETS[unlist(l1cfg$analysis_sets)]
MIN_OVERLAP <- as.integer(l1cfg$min_vmrs_overlapping)

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

## ----------------------------------------------------------- guard 2: model
repro <- rbindlist(lapply(regions, function(re) {
    r <- fit_hc3(feats[[re]], "line_l1_frac", covariates_for(re))
    s <- sealed[[re]]
    if (nrow(s) != 1L) stop("No unique sealed primary line_l1_frac row for ", re)
    data.table(region = re, estimate_refit = r$estimate, estimate_sealed = s$estimate,
               se_refit = r$se, se_sealed = s$se,
               n_fitted_refit = r$n_fitted, n_fitted_sealed = s$n_fitted)
}))
repro[, `:=`(abs_diff_estimate = abs(estimate_refit - estimate_sealed),
             abs_diff_se = abs(se_refit - se_sealed))]
print(repro)
if (any(!(repro$abs_diff_estimate <= REPRO_TOL & repro$abs_diff_se <= REPRO_TOL &
          repro$n_fitted_refit == repro$n_fitted_sealed))) {
    stop("The primary line_l1_frac fit does not reproduce the sealed cells.")
}
message("[04] guards passed: partition <= ", PARTITION_TOL, ", model <= ", REPRO_TOL)

## ------------------------------------------------------------------ the fits
assoc <- list(); contr <- list(); comp <- list()
for (re in regions) {
    d0 <- feats[[re]]
    for (sn in names(SETS)) {
        s <- SETS[[sn]]$subset(d0)
        ut <- usable_terms(s, c(covariates_for(re), SETS[[sn]]$extra))
        fits <- list()
        for (o in OUTCOMES) {
            n_any <- sum(s[[o]] > 0)
            comp[[length(comp) + 1]] <- data.table(
                region = re, analysis_set = sn, outcome = o, n_vmrs = nrow(s),
                n_vmrs_overlapping = n_any, mean_frac = mean(s[[o]]),
                mean_frac_among_overlapping = if (n_any) mean(s[[o]][s[[o]] > 0]) else NA_real_)
            if (n_any < MIN_OVERLAP) {
                assoc[[length(assoc) + 1]] <- data.table(
                    region = re, analysis_set = sn, outcome = o, n = nrow(s),
                    n_vmrs_overlapping = n_any, note = "insufficient_overlap")
                next
            }
            r <- fit_hc3(s, o, ut$keep)
            se_jk <- block_jackknife(
                s$chrom, function(keep) coef_pred(fit_glm(s[keep], o, ut$keep)), r$estimate)
            fits[[o]] <- r$estimate
            assoc[[length(assoc) + 1]] <- data.table(
                region = re, analysis_set = sn, outcome = o, n = nrow(s),
                n_vmrs_overlapping = n_any, n_fitted = r$n_fitted,
                estimate = r$estimate, se_hc3 = r$se, p_hc3 = r$p,
                se_jackknife = se_jk,
                p_jackknife = 2 * stats::pnorm(-abs(r$estimate / se_jk)),
                terms_dropped_in_set = paste(ut$dropped, collapse = ";"),
                note = NA_character_)
        }
        for (cn in names(CONTRASTS)) {
            a <- CONTRASTS[[cn]][[1]]; b <- CONTRASTS[[cn]][[2]]
            if (is.null(fits[[a]]) || is.null(fits[[b]])) {
                contr[[length(contr) + 1]] <- data.table(
                    region = re, analysis_set = sn, contrast = cn,
                    outcome_a = a, outcome_b = b, note = "insufficient_overlap")
                next
            }
            delta <- fits[[a]] - fits[[b]]
            se_joint <- block_jackknife(
                s$chrom,
                function(keep) {
                    kb <- s[keep]
                    coef_pred(fit_glm(kb, a, ut$keep)) - coef_pred(fit_glm(kb, b, ut$keep))
                },
                delta)
            contr[[length(contr) + 1]] <- data.table(
                region = re, analysis_set = sn, contrast = cn, outcome_a = a, outcome_b = b,
                estimate_a = fits[[a]], estimate_b = fits[[b]], difference = delta,
                se_joint_jackknife = se_joint, z = delta / se_joint,
                p = 2 * stats::pnorm(-abs(delta / se_joint)), note = NA_character_)
        }
    }
    message("[04] ", re, " subfamily fits done")
}
assoc <- rbindlist(assoc, fill = TRUE)
contr <- rbindlist(contr, fill = TRUE)
comp <- rbindlist(comp, fill = TRUE)

## ------------------------------------------------------- descriptive reading
rr <- l1cfg$reading_rule
reading <- rbindlist(lapply(names(CONTRASTS), function(cn) {
    pr <- contr[contrast == cn & analysis_set == "primary" &
                region %in% unlist(rr$claim_regions)]
    hm <- contr[contrast == cn & analysis_set == rr$require_sign_kept_in &
                region %in% unlist(rr$claim_regions)]
    ok_pr <- nrow(pr) == length(rr$claim_regions) && all(is.finite(pr$p)) &&
        all(pr$p < rr$alpha) && length(unique(sign(pr$difference))) == 1L
    s0 <- if (nrow(pr)) sign(pr$difference[1]) else NA_real_
    hm_fitted <- nrow(hm) == length(rr$claim_regions) && all(is.finite(hm$difference))
    ok_hm <- ok_pr && hm_fitted && all(sign(hm$difference) == s0)
    ## Three readings, so "could not be tested" is never reported as "tested and
    ## failed": a primary that passes but whose high-mappability arm falls under
    ## min_vmrs_overlapping is not_evaluable, not not_consistent.
    data.table(contrast = cn,
               primary_p_dlpfc = pr[region == "dlpfc", p][1],
               primary_p_hippocampus = pr[region == "hippocampus", p][1],
               primary_sign = if (ok_pr) s0 else NA_real_,
               high_mappability_fitted = hm_fitted,
               sign_kept_high_mappability = ok_hm,
               reading = if (ok_hm) "consistent" else if (ok_pr && !hm_fitted)
                   "not_evaluable_high_mappability" else "not_consistent")
}))

set_aside <- unlist(annot$interpretation$technically_confounded_regions$line_l1_frac)
label <- function(dt) {
    dt[, `:=`(set_aside_from_claim = region %in% set_aside,
              predictor = PRED, gating = FALSE, outside_bh_family = TRUE,
              sensitivity = "line_l1_subfamily_resolution",
              population = cohort,
              source_run_id = vapply(region, function(re) up[[re]]$run_id, character(1)),
              vmr_set_id = vapply(region, function(re) up[[re]]$vmr_set_id, character(1)))]
    dt[]
}
assoc <- label(assoc); contr <- label(contr); comp <- label(comp)

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
        config_l1_subfamilies_sha256 = attr(l1cfg, "config_sha256"),
        l1_subfamily_asset_sha256 = asset_sha,
        sensitivity = "line_l1_subfamily_resolution",
        gating = "FALSE",
        min_vmrs_overlapping = MIN_OVERLAP,
        reproduction_tolerance = REPRO_TOL,
        partition_tolerance = PARTITION_TOL
    )
)
res_dir <- file.path(run$dir, "results")
dir.create(res_dir, showWarnings = FALSE)
code_dir <- file.path(run$dir, "code", "_h")
dir.create(code_dir, recursive = TRUE)
invisible(file.copy(file.path(repo_root(), MODULE, "_h", "10_l1_subfamily.R"), code_dir))
cfg_dir <- file.path(run$dir, "code", "config")
dir.create(cfg_dir, recursive = TRUE)
invisible(file.copy(file.path(repo_root(), "config", c("l1_subfamilies.yml", "repeat_annotations.yml")),
                    cfg_dir))

write_atomic(asset_counts, file.path(res_dir, "l1-subfamily-asset-classes.tsv"))
write_atomic(partition, file.path(res_dir, "l1-subfamily-partition-check.tsv"))
write_atomic(repro, file.path(res_dir, "l1-subfamily-reproduction-check.tsv"))
write_atomic(comp, file.path(res_dir, "l1-subfamily-composition.tsv"))
write_atomic(assoc, file.path(res_dir, "l1-subfamily-association.tsv"))
write_atomic(contr, file.path(res_dir, "l1-subfamily-contrasts.tsv"))
write_atomic(reading, file.path(res_dir, "l1-subfamily-reading.tsv"))
write_atomic(feats_long <- rbindlist(lapply(regions, function(re)
    feats[[re]][, c("vmr_id", "chrom", "start", "end", "line_l1_frac", OUTCOMES), with = FALSE][
        , region := re])), file.path(res_dir, "l1-subfamily-vmr-features.tsv"))
writeLines(c(
    "Interpretation constraints carried by this run:",
    "  - NON-GATING and outside the BH family. No gate, claim, q-value or decision",
    "    token of the accepted cells changes.",
    "  - Overlap is overlap. Young or full-length L1 overlap does not show activity,",
    "    expression or retrotransposition, and full-length is not retrotransposition",
    "    competence (no ORF integrity check).",
    "  - Young L1 sequence is the least mappable in the genome. A young-L1 estimate",
    "    that loses its sign under high_mappability is not evidence.",
    "  - Full-length is conservative: elements split by an insertion are counted as",
    "    fragments (config/l1_subfamilies.yml:completeness).",
    "  - Caudate is fitted and shown but set aside from any claim, as caudate",
    "    line_l1_frac is (batch 3; AGENTS.md 8.1).",
    "  - Contrasts are differences of slopes on the same VMRs, tested with the JOINT",
    "    chromosome jackknife. 'consistent' in l1-subfamily-reading.tsv is a",
    "    descriptive reading rule fixed in config before fitting, not a gate."
), file.path(res_dir, "interpretation-constraints.txt"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "session-info.txt"))

if (!allow_unlocked) close_run(run)
message("[04] L1 subfamily run ", run$run_id, if (allow_unlocked) " (smoke, not sealed)" else " sealed")
print(contr[analysis_set %in% c("primary", "high_mappability"),
            .(region, analysis_set, contrast, estimate_a, estimate_b, difference,
              se_joint_jackknife, p)], digits = 3)
print(reading)
cat(run$run_id, "\n", sep = "")
