#!/usr/bin/env Rscript
#### 09 stage 17 -- is the SCZ axis depletion specific, or a property of GWAS loci? ####
##
## Usage:
##   Rscript _h/17_negative_control_traits.R --cohort AA [--allow-unlocked]
##   (after _h/17a_extract_gwas_leads.sh has written _m/combined/negative-control-gwas/sig/)
##
## Module 09 finds SCZ-linked VMRs LOWER on the local-genetic-control axis in
## every region, adjusted for technical covariates only. Module 04 shows the
## low-control end of the axis is gene-proximal and active; GWAS loci of most
## traits are gene-dense. This stage asks whether the depletion is specific to
## schizophrenia by running the identical locus -> VMR -> axis contrast for every
## trait in a harmonized GWAS collection under ONE locus rule, applied to
## schizophrenia too, and by refitting every contrast with broad genomic context
## as a covariate. Design and reading rules: config/gwas_negative_controls.yml.
##
## Two SCZ rows exist per region and must not be confused:
##   PGC3_SCZ            -- PGC3 under the uniform lead-SNP rule (comparable to
##                          every other trait; the row the distribution is read
##                          against);
##   PGC3_SCZ_published  -- the accepted Module 09 linkage (published fine-mapped
##                          intervals + flank), refitted here only so model B's
##                          context attenuation is reported on the claim's own
##                          linkage.
##
## Module-level: reads the accepted Module 09 runs and their upstreams; writes
## _m/combined/. It changes no decision.

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
## --gwas-dir: smoke use only, points the stage at a partial extraction.
gw_dir <- opts$gwas_dir %||% file.path(out_dir, "negative-control-gwas")
if (!file.exists(file.path(gw_dir, "sig", ".complete"))) {
    stop("Extraction incomplete: run _h/17a_extract_gwas_leads.sh first (", gw_dir, ")")
}
lr <- cfg$locus_rule; ax <- cfg$axis
predictor <- ax$predictor
tech <- as.character(unlist(ax$technical_covariates))
ctx <- as.character(ax$context_covariate)
flank <- as.integer(lr$link_flank_bp)
alpha <- as.numeric(ax$alpha)
confounded <- as.character(config_get(load_config("region_donor_generalization"),
    "interpretation.technically_confounded_regions.all_outcomes"))

## ------------------------------------------------------------- traits
meta <- fread(cfg$collection$metadata, colClasses = "character")
meta <- meta[, .(tag = Tag, category = Category, sample_size = as.numeric(Sample_Size),
                 binary = Binary == "1", keep = Keep)]
sig_files <- list.files(file.path(gw_dir, "sig"), pattern = "\\.tsv$", full.names = TRUE)
tags <- sub("\\.tsv$", "", basename(sig_files))
traits <- data.table(tag = tags, sig_file = sig_files)
traits <- merge(traits, meta, by = "tag", all.x = TRUE)

## Focal GWAS supplied outside the collection (config focal_sumstats): PGC3
## schizophrenia and the multi-ancestry releases. Their role, category and
## ancestry come from the config, so adding one needs no code change here.
focal <- rbindlist(lapply(cfg$focal_sumstats, function(e) data.table(
    tag = as.character(e$tag), focal_role = as.character(e$role),
    focal_category = as.character(e$category), ancestry = as.character(e$ancestry),
    focal_sample_size = as.numeric(e$sample_size %||% NA_real_),
    compare_to = as.character(e$compare_to %||% NA_character_))), fill = TRUE)
missing_focal <- setdiff(focal$tag, traits$tag)
if (length(missing_focal)) {
    stop("focal_sumstats not extracted (rerun _h/17a_extract_gwas_leads.sh): ",
         paste(missing_focal, collapse = ", "))
}
traits <- merge(traits, focal, by = "tag", all.x = TRUE)
traits[!is.na(focal_role), `:=`(category = focal_category, keep = "Yes", binary = TRUE,
                                sample_size = focal_sample_size)]
pos_ctl <- as.character(unlist(cfg$collection$positive_control_tags))
traits[, role := fcase(!is.na(focal_role), focal_role,
                       tag %in% pos_ctl, "positive_control",
                       isTRUE(cfg$collection$exclude_unless_keep_yes) & !keep %in% "Yes",
                           "excluded_not_keep",
                       default = "negative_control")]

## ------------------------------------------------------------- lead clumping
excl <- rbindlist(lapply(lr$exclude_regions_hg38, as.data.table))
clump_leads <- function(f) {
    s <- fread(f, colClasses = list(character = c("chrom", "variant_id")))
    s[, chrom := paste0("chr", sub("^chr", "", chrom))]
    s <- s[is.finite(pvalue) & is.finite(pos)]
    if (isTRUE(lr$autosomes_only)) s <- s[chrom %in% paste0("chr", 1:22)]
    for (i in seq_len(nrow(excl))) {
        s <- s[!(chrom == excl$chrom[i] & pos >= excl$start[i] & pos <= excl$end[i])]
    }
    n_sig <- nrow(s)
    if (!n_sig) return(list(leads = s[0], n_sig = 0L))
    setorder(s, pvalue)
    d <- as.integer(lr$clump_distance_bp)
    leads <- list()
    for (ch in unique(s$chrom)) {
        sc <- s[chrom == ch]
        kept <- integer(0)
        for (i in seq_len(nrow(sc))) {
            if (!length(kept) || all(abs(sc$pos[i] - sc$pos[kept]) > d)) kept <- c(kept, i)
        }
        leads[[ch]] <- sc[kept]
    }
    list(leads = rbindlist(leads), n_sig = n_sig)
}

## ------------------------------------------------------------- per region
vmr_tab <- list(); feat_tab <- list(); scz_link <- list(); run_ids <- list()
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
    need <- c(tech, ctx)
    miss <- setdiff(need, names(feat))
    if (length(miss)) stop("Module 04 features lack: ", paste(miss, collapse = ","))
    v <- merge(v, feat[, c("vmr_id", need), with = FALSE], by = "vmr_id", all.x = TRUE)
    v[[ctx]] <- factor(v[[ctx]])
    vmr_tab[[re]] <- v
    sl <- fread(file.path(module_root, "_m", "runs", r09$run_id, "results", "vmr-scz-linkage.tsv"))
    scz_link[[re]] <- sl[, .(vmr_id, linked = as.logical(scz_linked))]
}

fit_contrast <- function(v, linked) {
    d <- copy(v); d[, linked := linked]
    d <- d[is.finite(get(predictor)) & complete.cases(d[, c(tech, ctx), with = FALSE])]
    one <- function(covs, label) {
        if (length(unique(d$linked)) < 2) {
            return(data.table(model = label, estimate = NA_real_, se = NA_real_, p = NA_real_))
        }
        form <- stats::reformulate(c(predictor, covs), response = "linked")
        fit <- suppressWarnings(stats::glm(form, data = d, family = stats::binomial()))
        co <- summary(fit)$coefficients
        data.table(model = label, estimate = co[predictor, "Estimate"],
                   se = co[predictor, "Std. Error"], p = co[predictor, "Pr(>|z|)"])
    }
    rbind(one(tech, "A_technical"), one(c(tech, ctx), "B_technical_context"))[,
        `:=`(n_linked = sum(d$linked), n_background = sum(!d$linked))]
}

rows <- list(); lead_rows <- list()
for (i in seq_len(nrow(traits))) {
    tr <- traits[i]
    cl <- clump_leads(tr$sig_file)
    leads <- cl$leads
    lead_rows[[tr$tag]] <- if (nrow(leads)) leads[, .(tag = tr$tag, chrom, pos, pvalue, variant_id)] else NULL
    for (re in regions) {
        v <- vmr_tab[[re]]
        linked <- rep(FALSE, nrow(v)); n_leads_with_vmr <- 0L
        if (nrow(leads)) {
            w <- leads[, .(chrom, start = pmax(1L, pos - flank), end = pos + flank, lead = .I)]
            setkey(w, chrom, start, end)
            ov <- foverlaps(v[, .(vmr_id, chrom, start, end)], w,
                            by.x = c("chrom", "start", "end"), type = "any", nomatch = NULL)
            linked <- v$vmr_id %in% ov$vmr_id
            n_leads_with_vmr <- uniqueN(ov$lead)
        }
        res <- fit_contrast(v, linked)
        res[, `:=`(tag = tr$tag, role = tr$role, category = tr$category,
                   sample_size = tr$sample_size, binary = tr$binary,
                   ancestry = tr$ancestry, region = re,
                   locus_rule = "uniform_lead_snp", n_sig_variants = cl$n_sig,
                   n_leads = nrow(leads), n_leads_with_vmr = n_leads_with_vmr,
                   n_vmrs = nrow(v))]
        rows[[length(rows) + 1]] <- res
    }
    if (i %% 10 == 0) message("[17] ", i, "/", nrow(traits), " traits")
}
## The claim's own linkage (published PGC3 intervals), for the context arm.
for (re in regions) {
    v <- vmr_tab[[re]]
    lk <- scz_link[[re]]$linked[match(v$vmr_id, scz_link[[re]]$vmr_id)]
    lk[is.na(lk)] <- FALSE
    res <- fit_contrast(v, lk)
    res[, `:=`(tag = "PGC3_SCZ_published", role = "scz_published_intervals",
               category = "Psychiatric-neurologic", sample_size = NA_real_, binary = TRUE,
               ancestry = "EUR",
               region = re, locus_rule = "pgc3_fine_mapped_intervals_plus_flank",
               n_sig_variants = NA_integer_, n_leads = NA_integer_,
               n_leads_with_vmr = NA_integer_, n_vmrs = nrow(v))]
    rows[[length(rows) + 1]] <- res
}
all <- rbindlist(rows, fill = TRUE)
all[, frac_vmrs_linked := n_linked / (n_linked + n_background)]
all[, in_distribution := role == "negative_control" & n_leads >= as.integer(lr$min_leads_for_distribution)]
all[, q := NA_real_]
all[in_distribution == TRUE, q := stats::p.adjust(p, ax$fdr_method), by = .(region, model)]
all[, `:=`(cohort = opts$cohort, tier = fifelse(region %in% confounded, "descriptive_only", "claim_eligible"),
           module_09_run_id = unlist(run_ids)[region])]

## ------------------------------------------------------------- reading
scz_u <- all[tag == "PGC3_SCZ"]
summ <- rbindlist(lapply(regions, function(re) rbindlist(lapply(c("A_technical", "B_technical_context"), function(m) {
    s <- scz_u[region == re & model == m]
    d <- all[region == re & model == m & in_distribution == TRUE & is.finite(estimate)]
    sim <- d[n_leads >= s$n_leads / as.numeric(ax$similar_lead_count_factor) &
             n_leads <= s$n_leads * as.numeric(ax$similar_lead_count_factor)]
    psy <- d[category == "Psychiatric-neurologic"]
    nonpsy <- d[category != "Psychiatric-neurologic"]
    pub <- all[region == re & model == m & tag == "PGC3_SCZ_published"]
    data.table(
        region = re, model = m,
        scz_estimate = s$estimate, scz_se = s$se, scz_p = s$p, scz_n_leads = s$n_leads,
        scz_published_estimate = pub$estimate, scz_published_p = pub$p,
        n_traits_in_distribution = nrow(d),
        distribution_median = stats::median(d$estimate),
        distribution_q25 = stats::quantile(d$estimate, 0.25),
        distribution_q75 = stats::quantile(d$estimate, 0.75),
        n_traits_negative_q_below_alpha = sum(d$q < alpha & d$estimate < 0),
        n_traits_positive_q_below_alpha = sum(d$q < alpha & d$estimate > 0),
        scz_percentile_in_distribution = mean(d$estimate < s$estimate),
        n_traits_more_negative_than_scz = sum(d$estimate < s$estimate),
        similar_lead_count_n = nrow(sim),
        similar_lead_count_median = if (nrow(sim)) stats::median(sim$estimate) else NA_real_,
        scz_percentile_similar_lead_count = if (nrow(sim)) mean(sim$estimate < s$estimate) else NA_real_,
        psychiatric_median = if (nrow(psy)) stats::median(psy$estimate) else NA_real_,
        nonpsychiatric_median = if (nrow(nonpsy)) stats::median(nonpsy$estimate) else NA_real_,
        psychiatric_vs_other_wilcoxon_p = if (nrow(psy) > 1 && nrow(nonpsy) > 1)
            suppressWarnings(stats::wilcox.test(psy$estimate, nonpsy$estimate)$p.value) else NA_real_,
        spearman_estimate_vs_log_n_leads = suppressWarnings(stats::cor(d$estimate, log(d$n_leads), method = "spearman")),
        spearman_estimate_vs_frac_linked = suppressWarnings(stats::cor(d$estimate, d$frac_vmrs_linked, method = "spearman")))
}))))
summ[, scz_context_attenuation := NA_real_]
for (re in regions) {
    a <- summ[region == re & model == "A_technical", scz_published_estimate]
    b <- summ[region == re & model == "B_technical_context", scz_published_estimate]
    summ[region == re, scz_context_attenuation := 1 - b / a]
}
summ[, `:=`(cohort = opts$cohort, tier = fifelse(region %in% confounded, "descriptive_only", "claim_eligible"),
            post_hoc_sensitivity = TRUE, alters_module_09_decision = FALSE,
            regions_are_independent_replicates = FALSE,
            built_with_unaccepted_runs = allow_unaccepted)]
## ------------------------------------------------- ancestry sensitivity
## The GWAS collection is European-dominated and the methylation cohort is
## admixed African American. No African-ancestry psychiatric GWAS yields a single
## genome-wide significant locus (config ancestry:), so a matched test is
## impossible and these multi-ancestry releases are the available partial check.
## Each is placed against the null distribution and, where the config names one,
## against its own European counterpart under the identical rule.
anc <- rbindlist(lapply(regions, function(re) rbindlist(lapply(
    c("A_technical", "B_technical_context"), function(m) {
    d <- all[region == re & model == m & in_distribution == TRUE & is.finite(estimate)]
    rbindlist(lapply(traits[role == "ancestry_sensitivity"]$tag, function(tg) {
        s <- all[region == re & model == m & tag == tg]
        cmp_tag <- traits[tag == tg]$compare_to[1]
        cmp <- if (!is.na(cmp_tag)) all[region == re & model == m & tag == cmp_tag] else NULL
        data.table(
            region = re, model = m, tag = tg,
            ancestry = traits[tag == tg]$ancestry[1],
            estimate = s$estimate, se = s$se, p = s$p, n_leads = s$n_leads,
            percentile_in_distribution = if (nrow(d)) mean(d$estimate < s$estimate) else NA_real_,
            distribution_median = if (nrow(d)) stats::median(d$estimate) else NA_real_,
            compared_with = cmp_tag %||% NA_character_,
            comparator_estimate = if (!is.null(cmp) && nrow(cmp)) cmp$estimate else NA_real_,
            comparator_n_leads = if (!is.null(cmp) && nrow(cmp)) cmp$n_leads else NA_integer_,
            ## Positive = the multi-ancestry release gives a WEAKER depletion.
            shift_from_comparator = if (!is.null(cmp) && nrow(cmp)) s$estimate - cmp$estimate else NA_real_)
    }), fill = TRUE)
}))))
anc[, `:=`(cohort = opts$cohort,
           tier = fifelse(region %in% confounded, "descriptive_only", "claim_eligible"),
           ## The African-ancestry samples are a small minority of each release,
           ## so this is a partial check and never an ancestry-matched one.
           matched_ancestry_test = FALSE,
           built_with_unaccepted_runs = allow_unaccepted)]

by_cat <- all[in_distribution == TRUE & is.finite(estimate),
              .(n_traits = .N, median_estimate = stats::median(estimate),
                n_negative_q_below_alpha = sum(q < alpha & estimate < 0),
                n_positive_q_below_alpha = sum(q < alpha & estimate > 0)),
              by = .(region, model, category)][order(region, model, median_estimate)]

sfx <- if (allow_unaccepted) paste0("-", opts$cohort, "-UNACCEPTED") else paste0("-", opts$cohort)
write_atomic(all, file.path(out_dir, paste0("scz-negative-control-traits", sfx, ".tsv")))
write_atomic(summ, file.path(out_dir, paste0("scz-negative-control-summary", sfx, ".tsv")))
write_atomic(by_cat, file.path(out_dir, paste0("scz-negative-control-by-category", sfx, ".tsv")))
write_atomic(anc, file.path(out_dir, paste0("scz-negative-control-ancestry", sfx, ".tsv")))
write_atomic(rbindlist(lead_rows), file.path(out_dir, paste0("scz-negative-control-leads", sfx, ".tsv")))
options(width = 220)
print(summ[, .(region, model, scz_estimate, scz_p, scz_n_leads, n_traits_in_distribution,
               distribution_median, scz_percentile_in_distribution,
               n_traits_negative_q_below_alpha, scz_percentile_similar_lead_count,
               psychiatric_median, nonpsychiatric_median, scz_context_attenuation)])
print(anc[, .(region, model, tag, n_leads, estimate, p, percentile_in_distribution,
              compared_with, comparator_estimate, shift_from_comparator)])
message("[17] negative-control traits written to ", out_dir,
        if (allow_unaccepted) "  (UNACCEPTED: unlocked config)" else "")
