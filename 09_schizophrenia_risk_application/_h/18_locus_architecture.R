#!/usr/bin/env Rscript
#### 09 stage 18 -- what ARE the trait loci, in Module 04 architecture terms? ####
##
## Usage:
##   Rscript _h/18_locus_architecture.R --cohort AA [--allow-unlocked]
##   (after _h/17_negative_control_traits.R)
##
## Stage 17 establishes that schizophrenia is a typical GWAS trait on the
## local-genetic-control axis: most well-powered traits show the same depletion.
## It does not say what separates the traits that ARE depleted from the ones that
## are not, and that separation is the biological question. Two facts from the
## 2026-09-19 stage 17 run motivate this stage, both recorded in
## config/gwas_negative_controls.yml before it was written:
##   - anthropometric traits are the most depleted category in all three regions;
##   - immune traits show no depletion at all and trend positive outside caudate.
##
## So: do immune loci sit in repeat-rich, heterochromatic, high-control sequence
## while anthropometric loci sit in gene-dense active sequence? If yes, the axis
## contrast is reporting WHERE a trait's loci live in the genome, which ties the
## negative control to the Module 04 architecture (AGENTS.md 7.4) instead of
## merely qualifying Module 09.
##
## Design (config locus_architecture), mirroring 09b / AGENTS.md 7.9:
##   - ONE annotation at a time, never a joint annotation model;
##   - each fitted twice, with and without the control score, so an enrichment
##     that is only the axis restated is visible as such;
##   - a high-mappability arm for the repeat and heterochromatin annotations,
##     which AGENTS.md 7.4 requires before any such statement;
##   - same glm, same technical covariates and same linkage as stage 17, so the
##     estimates sit on one scale with the axis estimate.
##
## Estimate = log odds ratio for a VMR near a trait lead carrying the annotation.
## DESCRIPTIVE. Overlap is overlap: no activity, retrotransposition, cell-type or
## trait-mechanism claim follows from it.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

`%||%` <- function(a, b) if (is.null(a)) b else a
MODULE <- "09_schizophrenia_risk_application"
opts <- parse_v2_args(require = "cohort")
allow_unaccepted <- isTRUE(opts$allow_unlocked)
cfg <- load_config("gwas_negative_controls")
assert_locked(list(gwas_negative_controls = cfg), allow_unlocked = allow_unaccepted)
scz <- load_config("schizophrenia")
regions <- as.character(unlist(scz$regions %||% c("caudate", "dlpfc", "hippocampus")))
module_root <- file.path(repo_root(), MODULE)
out_dir <- file.path(module_root, "_m", "combined")

lr <- cfg$locus_rule; ax <- cfg$axis; la <- cfg$locus_architecture
predictor <- ax$predictor
tech <- as.character(unlist(ax$technical_covariates))
flank <- as.integer(lr$link_flank_bp)
alpha <- as.numeric(ax$alpha)
indicators <- as.character(unlist(la$indicators))
continuous <- as.character(unlist(la$continuous_z))
annotations <- c(indicators, continuous)
min_linked <- as.integer(la$min_linked_vmrs)
mr <- la$mappability_restricted_arm
mr_anns <- as.character(unlist(mr$annotations))
mr_min <- as.numeric(mr$min_mappability)
confounded <- as.character(config_get(load_config("region_donor_generalization"),
    "interpretation.technically_confounded_regions.all_outcomes"))

sfx <- if (allow_unaccepted) paste0("-", opts$cohort, "-UNACCEPTED") else paste0("-", opts$cohort)
## Stage 17's own outputs are the input: its clumped leads and its axis estimates.
leads_file <- file.path(out_dir, paste0("scz-negative-control-leads", sfx, ".tsv"))
traits_file <- file.path(out_dir, paste0("scz-negative-control-traits", sfx, ".tsv"))
for (f in c(leads_file, traits_file)) {
    if (!file.exists(f)) stop("Stage 17 output missing, run _h/17_negative_control_traits.R first: ", f)
}
leads_all <- fread(leads_file, colClasses = list(character = c("chrom", "variant_id")))
axis_tab <- fread(traits_file)
## One row per trait: role, category and the stage 17 axis estimate to read
## annotation enrichment against. Model A, the technical-covariate fit.
trait_meta <- unique(axis_tab[, .(tag, role, category, ancestry, n_leads, in_distribution)])
axis_a <- axis_tab[model == "A_technical", .(tag, region, axis_estimate = estimate,
                                             axis_p = p, n_linked_axis = n_linked)]

## ------------------------------------------------------------- per region
vmr_tab <- list(); run_ids <- list()
for (re in regions) {
    r09 <- require_accepted_upstream(MODULE, opts$cohort, re, allow_unaccepted = allow_unaccepted)
    run_ids[[re]] <- r09$run_id
    man <- fread(file.path(module_root, "_m", "runs", r09$run_id, "manifest.tsv"),
                 colClasses = "character")
    mf <- function(field) man$value[man$field == field][1]
    lgc <- load_local_genetic_control(mf("upstream_local_genetic_variance_run_id"),
                                      region = re, cohort = opts$cohort, eligible_only = TRUE)
    v <- as.data.table(lgc)[, c("vmr_id", "chrom", "start", "end", predictor), with = FALSE]
    v[, chrom := paste0("chr", sub("^chr", "", chrom))]
    feat <- fread(file.path(repo_root(), "04_repeat_repressive_architecture", "_m", "runs",
                            mf("upstream_repeat_architecture_run_id"), "results",
                            "vmr-features.tsv"))
    banned <- intersect(unlist(scz$forbidden_columns), names(feat))
    if (length(banned)) stop("Module 04 features carry forbidden column(s): ",
                             paste(banned, collapse = ","))
    need <- unique(c(tech, annotations))
    miss <- setdiff(need, names(feat))
    if (length(miss)) stop("Module 04 features lack: ", paste(miss, collapse = ","))
    v <- merge(v, feat[, c("vmr_id", need), with = FALSE], by = "vmr_id", all.x = TRUE)
    ## Indicators to 0/1; continuous annotations standardized WITHIN region, so a
    ## coefficient is per SD and regions are never compared on a raw scale.
    for (cc in indicators) v[[cc]] <- as.numeric(as.logical(v[[cc]]))
    for (cc in continuous) {
        x <- as.numeric(v[[cc]])
        v[[cc]] <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
    }
    vmr_tab[[re]] <- v
    message("[18] ", re, ": ", nrow(v), " eligible VMRs, run ", r09$run_id)
}

## ------------------------------------------------------------- fitting
## glm.fit on a prebuilt model matrix: identical coefficients and SEs to glm(),
## without re-parsing a formula ~10,000 times. Asserted against glm() below.
fit_one <- function(X, y) {
    ok <- stats::complete.cases(X) & !is.na(y)
    if (sum(ok) < 20L || length(unique(y[ok])) < 2L || length(unique(X[ok, 2])) < 2L) {
        return(c(NA_real_, NA_real_, NA_real_))
    }
    f <- try(suppressWarnings(stats::glm.fit(X[ok, , drop = FALSE], y[ok],
                                             family = stats::binomial())), silent = TRUE)
    if (inherits(f, "try-error") || !isTRUE(f$converged)) return(c(NA_real_, NA_real_, NA_real_))
    est <- f$coefficients[2]
    se <- sqrt(diag(chol2inv(qr.R(f$qr))))[2]
    if (!is.finite(est) || !is.finite(se) || se <= 0) return(c(NA_real_, NA_real_, NA_real_))
    c(est, se, 2 * stats::pnorm(-abs(est / se)))
}

## Self-check: the fast path must reproduce glm() on a real design.
local({
    v <- vmr_tab[[regions[1]]]
    y <- as.numeric(v[[indicators[1]]])
    X <- cbind(1, v[[predictor]], as.matrix(v[, tech, with = FALSE]))
    a <- fit_one(X, y)
    d <- data.frame(y = y, p = v[[predictor]], v[, tech, with = FALSE])
    b <- summary(suppressWarnings(stats::glm(y ~ ., data = d,
        family = stats::binomial())))$coefficients[2, 1:2]
    if (!isTRUE(all.equal(unname(a[1:2]), unname(b), tolerance = 1e-8))) {
        stop("glm.fit fast path disagrees with glm(): ", paste(signif(a[1:2], 10), collapse = ","),
             " vs ", paste(signif(b, 10), collapse = ","))
    }
    message("[18] fast-path check passed against glm()")
})

## Arm definitions: which VMR rows each arm uses.
arm_rows <- function(v, arm) {
    if (arm == "all_vmrs") rep(TRUE, nrow(v))
    else !is.na(v$mappability) & v$mappability >= mr_min
}
arms_for <- function(ann) {
    if (ann %in% mr_anns) c("all_vmrs", "high_mappability") else "all_vmrs"
}
adjustments <- names(la$adjustments)

rows <- list()
tags <- sort(unique(leads_all$tag))
for (re in regions) {
    v <- vmr_tab[[re]]
    vk <- v[, .(vmr_id, chrom, start, end)]
    ## Prebuild, per arm x adjustment, the covariate block that does not change
    ## across traits or annotations.
    base <- list()
    for (arm in c("all_vmrs", "high_mappability")) {
        keep <- arm_rows(v, arm)
        for (adj in adjustments) {
            extra <- as.character(unlist(la$adjustments[[adj]]))
            cols <- c(tech, extra)
            base[[paste(arm, adj)]] <- list(
                keep = keep,
                C = as.matrix(v[keep, cols, with = FALSE]))
        }
    }
    for (tg in tags) {
        ld <- leads_all[tag == tg]
        if (!nrow(ld)) next
        w <- ld[, .(chrom, start = pmax(1L, pos - flank), end = pos + flank)]
        setkey(w, chrom, start, end)
        ov <- foverlaps(vk, w, by.x = c("chrom", "start", "end"), type = "any", nomatch = NULL)
        linked <- as.numeric(v$vmr_id %in% ov$vmr_id)
        if (sum(linked) < min_linked) {
            rows[[length(rows) + 1L]] <- data.table(
                region = re, tag = tg, annotation = NA_character_, adjustment = NA_character_,
                arm = NA_character_, estimate = NA_real_, se = NA_real_, p = NA_real_,
                n_linked = sum(linked), n_background = sum(!linked),
                mean_linked = NA_real_, mean_background = NA_real_,
                skipped = TRUE, skip_reason = paste0("fewer_than_", min_linked, "_linked_vmrs"))
            next
        }
        for (ann in annotations) {
            a <- as.numeric(v[[ann]])
            for (arm in arms_for(ann)) {
                for (adj in adjustments) {
                    b <- base[[paste(arm, adj)]]
                    keep <- b$keep
                    X <- cbind(1, a[keep], b$C)
                    r <- fit_one(X, linked[keep])
                    rows[[length(rows) + 1L]] <- data.table(
                        region = re, tag = tg, annotation = ann, adjustment = adj, arm = arm,
                        estimate = r[1], se = r[2], p = r[3],
                        n_linked = sum(linked[keep]), n_background = sum(!linked[keep]),
                        mean_linked = mean(a[keep][linked[keep] == 1], na.rm = TRUE),
                        mean_background = mean(a[keep][linked[keep] == 0], na.rm = TRUE),
                        skipped = FALSE, skip_reason = NA_character_)
                }
            }
        }
    }
    message("[18] ", re, ": fitted ", length(tags), " traits")
}
arch <- rbindlist(rows, fill = TRUE)
arch <- merge(arch, trait_meta, by = "tag", all.x = TRUE)
arch <- merge(arch, axis_a, by = c("tag", "region"), all.x = TRUE)
arch[, q := NA_real_]
arch[in_distribution == TRUE & skipped == FALSE & is.finite(p),
     q := stats::p.adjust(p, la$fdr_method), by = .(region, annotation, adjustment, arm)]
arch[, `:=`(cohort = opts$cohort,
            tier = fifelse(region %in% confounded, "descriptive_only", "claim_eligible"),
            post_hoc_sensitivity = TRUE,
            module_09_run_id = unlist(run_ids)[region],
            ## AGENTS.md 7.4: overlap licenses none of these.
            retrotransposition_claim_allowed = FALSE,
            cell_type_claim_allowed = FALSE,
            trait_mechanism_claim_allowed = FALSE,
            built_with_unaccepted_runs = allow_unaccepted)]

## ------------------------------------------------- does architecture track the axis?
## Over traits in the null distribution: a trait whose loci carry more of an
## annotation -- does it sit lower or higher on the control axis? This is the
## quantity that would connect stage 17 to the Module 04 architecture.
axis_link <- arch[in_distribution == TRUE & skipped == FALSE &
                  is.finite(estimate) & is.finite(axis_estimate),
                  {
                      n <- .N
                      rho <- if (n >= 10) suppressWarnings(stats::cor(estimate, axis_estimate,
                                             method = "spearman")) else NA_real_
                      pv <- if (n >= 10) suppressWarnings(stats::cor.test(estimate, axis_estimate,
                                             method = "spearman")$p.value) else NA_real_
                      .(n_traits = n, spearman_rho = rho, spearman_p = pv,
                        median_annotation_estimate = stats::median(estimate))
                  },
                  by = .(region, annotation, adjustment, arm)]
axis_link[, q := NA_real_]
axis_link[is.finite(spearman_p), q := stats::p.adjust(spearman_p, la$fdr_method),
          by = .(region, adjustment, arm)]
axis_link[, `:=`(cohort = opts$cohort,
                 tier = fifelse(region %in% confounded, "descriptive_only", "claim_eligible"),
                 interpretation = "negative_rho_means_annotation_marks_axis_depleted_traits")]

## ------------------------------------------------- per category
by_cat <- arch[in_distribution == TRUE & skipped == FALSE & is.finite(estimate),
               .(n_traits = .N, median_estimate = stats::median(estimate),
                 n_positive_q_below_alpha = sum(q < alpha & estimate > 0, na.rm = TRUE),
                 n_negative_q_below_alpha = sum(q < alpha & estimate < 0, na.rm = TRUE)),
               by = .(region, annotation, adjustment, arm, category)]
by_cat[, `:=`(cohort = opts$cohort,
              tier = fifelse(region %in% confounded, "descriptive_only", "claim_eligible"))]

write_atomic(arch, file.path(out_dir, paste0("scz-locus-architecture", sfx, ".tsv")))
write_atomic(axis_link, file.path(out_dir, paste0("scz-locus-architecture-axis-link", sfx, ".tsv")))
write_atomic(by_cat, file.path(out_dir, paste0("scz-locus-architecture-by-category", sfx, ".tsv")))

options(width = 220)
message("\n[18] annotation vs axis estimate across traits (score-adjusted, all VMRs):")
print(axis_link[adjustment == "technical_score" & arm == "all_vmrs"][order(region, spearman_rho)])
message("\n[18] category medians, score-adjusted, all VMRs, non-caudate regions:")
print(dcast(by_cat[adjustment == "technical_score" & arm == "all_vmrs" & region != "caudate" &
                   category %in% c("Immune", "Anthropometric", "Psychiatric-neurologic", "Blood")],
            annotation + category ~ region, value.var = "median_estimate"), digits = 2)
message("\n[18] schizophrenia and the multi-ancestry releases (score-adjusted, all VMRs):")
print(arch[tag %in% c("PGC3_SCZ", "PGC3_SCZ_multiancestry", "BIP2024_multiancestry") &
           adjustment == "technical_score" & arm == "all_vmrs" & skipped == FALSE,
           .(region, tag, annotation, estimate, p, mean_linked, mean_background)][order(region, tag, annotation)])
message("[18] locus architecture written to ", out_dir,
        if (allow_unaccepted) "  (UNACCEPTED: unlocked config)" else "")
