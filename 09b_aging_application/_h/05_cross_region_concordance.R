#!/usr/bin/env Rscript
#### 09b_aging_application stage 05 -- cross-region reading ####
##
## Usage:
##   Rscript _h/05_cross_region_concordance.R --cohort AA
##   Rscript _h/05_cross_region_concordance.R --cohort AA --allow-unlocked \
##       [--runs caudate=age-...,dlpfc=age-...,hippocampus=age-...]
##
## Every other stage runs inside ONE region. This one answers the cross-region
## half of the question, in the Module 08 tiers, and resolves the module's
## decision token. It asks three things and never pools them:
##
##   Q1 REGION-GENERAL (tier 1 logic, Module 09's H1 rule): is the association
##      supported in >= 2 regions, at least one of them not caudate? Endpoint is
##      concordance of per-region readings. No pooled p: the regions share
##      donors (DLPFC/hippocampus 115 of 118), so they are tissues within one
##      cohort, not independent replicates.
##
##   Q2 IDENTIFIED DIFFERENCE (tier 2): does the axis coefficient differ between
##      DLPFC and hippocampus, the one region contrast inside one sequencing
##      batch? Variance = PAIRED donor bootstrap + JOINT chromosome jackknife.
##      An independent-SE z would be wrong here: with nearly identical donors
##      the two coefficients are strongly positively correlated, and the paired
##      bootstrap measures that correlation instead of assuming it away. The
##      relative outcome makes the coefficients unit-free and commensurable; the
##      raw score is never compared across regions (AGENTS.md 7.2).
##
##   Q3 CAUDATE (tier 4): fitted and surfaced, excluded from any contrast.
##      Caudate is sequencing batch 3 (AGENTS.md 8.1).
##
## It also collates the DESCRIPTIVE annotation associations (stage 02) into one
## table and states, per annotation, in how many regions it is associated with
## age at q < alpha and whether those regions agree in sign. The same tier rule
## applies: consistency is counted outside caudate, caudate is shown beside it.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
H_DIR <- Sys.getenv("V2_RUN_CODE",
                    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09b_aging_application", "_h"))
source(file.path(H_DIR, "age_functions.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "09b_aging_application"
opts <- parse_v2_args(require = "cohort")
allow_unaccepted <- isTRUE(opts$allow_unlocked)
cfg <- load_config("aging")
xr <- cfg$cross_region
module_root <- file.path(repo_root(), MODULE)
out_dir <- file.path(module_root, "_m", "combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
regions <- as.character(unlist(cfg$regions))
confounded <- as.character(config_get(load_config("region_donor_generalization"),
    "interpretation.technically_confounded_regions.all_outcomes"))

## ------------------------------------------------------------- the runs
override <- list()
if (!is.null(opts$runs)) {
    if (!allow_unaccepted) stop("--runs is for smoke use only (--allow-unlocked)")
    for (kv in strsplit(opts$runs, ",", fixed = TRUE)[[1]]) {
        p <- strsplit(kv, "=", fixed = TRUE)[[1]]
        override[[p[1]]] <- p[2]
    }
}
resolve_run <- function(region) {
    if (!is.null(override[[region]])) return(override[[region]])
    if (!allow_unaccepted) {
        return(require_accepted_upstream(MODULE, opts$cohort, region)$run_id)
    }
    acc <- tryCatch(require_accepted_upstream(MODULE, opts$cohort, region)$run_id,
                    error = function(e) NA_character_)
    if (!is.na(acc)) return(acc)
    cand <- list.files(file.path(module_root, "_m", "runs"),
                       pattern = paste0("^age-(smoke-)?", opts$cohort, "-", region, "-"))
    if (!length(cand)) stop("No 09b run found for ", region)
    sort(cand, decreasing = TRUE)[[1L]]
}
runs <- vapply(regions, resolve_run, character(1))
run_dir <- function(re) file.path(module_root, "_m", "runs", runs[[re]])

read_run <- function(re) {
    dec <- fread(file.path(run_dir(re), "results", "aging-decision.tsv"))
    man <- fread(file.path(run_dir(re), "manifest.tsv"), colClasses = "character")
    cited <- man$value[man$field == "upstream_local_genetic_variance_run_id"]
    acc02 <- tryCatch(
        require_accepted_upstream("02_local_genetic_variance", opts$cohort, re)$run_id,
        error = function(e) NA_character_)
    dec[, `:=`(module_02_cited = cited, module_02_accepted = acc02,
               upstream_current = identical(cited, acc02),
               tier = if (re %in% confounded) "descriptive_only" else "claim_eligible")]
    dec
}
per_region <- rbindlist(lapply(regions, read_run), fill = TRUE)

## ------------------------------------------------------------- annotation collation
ann_all <- rbindlist(lapply(regions, function(re) {
    f <- file.path(run_dir(re), "results", "annotation-age-associations.tsv")
    if (!file.exists(f)) return(NULL)
    fread(f)
}), fill = TRUE)
ann_summary <- NULL
if (nrow(ann_all)) {
    a_alpha <- as.numeric(cfg$inference$alpha)
    tested <- ann_all[is.finite(estimate)]
    tested[, tier := fifelse(region %in% confounded, "descriptive_only", "claim_eligible")]
    ann_summary <- tested[, {
        hit <- q < a_alpha
        nc <- tier == "claim_eligible"
        signs <- unique(sign(estimate[hit & nc]))
        .(n_regions_tested = .N,
          regions_q_below_alpha = paste(region[hit], collapse = ","),
          n_noncaudate_q_below_alpha = sum(hit & nc),
          noncaudate_signs_agree = if (sum(hit & nc) >= 2L) length(signs) == 1L else NA,
          caudate_q = if (any(!nc)) q[!nc][1] else NA_real_,
          caudate_estimate = if (any(!nc)) estimate[!nc][1] else NA_real_,
          dlpfc_estimate = if (any(region == "dlpfc")) estimate[region == "dlpfc"][1] else NA_real_,
          dlpfc_q = if (any(region == "dlpfc")) q[region == "dlpfc"][1] else NA_real_,
          hippocampus_estimate = if (any(region == "hippocampus")) estimate[region == "hippocampus"][1] else NA_real_,
          hippocampus_q = if (any(region == "hippocampus")) q[region == "hippocampus"][1] else NA_real_)
    }, by = .(annotation, source_module, adjustment, outcome)]
    ann_summary[, `:=`(role = "descriptive", regions_are_independent_replicates = FALSE,
                       pooled_p_emitted = FALSE)]
}

## ------------------------------------------------------------- donor overlap
ck <- lapply(regions, function(re) readRDS(file.path(run_dir(re), "checkpoint",
                                                     "age-inputs.rds")))
names(ck) <- regions
donors <- lapply(ck, function(k) as.character(k$donors$FID))
overlap <- rbindlist(lapply(utils::combn(regions, 2L, simplify = FALSE), function(p) {
    a <- donors[[p[1]]]; b <- donors[[p[2]]]
    data.table(region_a = p[1], region_b = p[2], n_a = length(a), n_b = length(b),
               n_shared = length(intersect(a, b)),
               jaccard = length(intersect(a, b)) / length(union(a, b)))
}))

## ------------------------------------------------------------- Q1 region-general
supported <- per_region$region_supported %in% TRUE
outside <- !per_region$region %in% as.character(unlist(xr$requires_support_outside))
n_supported <- sum(supported)
region_general <- n_supported >= as.integer(xr$min_regions_supported) &&
    sum(supported & outside) >= 1L

## ------------------------------------------------------------- Q2 identified difference
contrast <- as.character(unlist(xr$identified_difference_contrast))
if (length(intersect(contrast, confounded))) {
    stop("The identified-difference contrast includes a technically confounded ",
         "region; tier 2 exists because it avoids the batch boundary.")
}
B <- if (allow_unaccepted) as.integer(xr$smoke_bootstrap_n) else as.integer(xr$bootstrap_n)
z_ci <- stats::qnorm(1 - (1 - as.numeric(xr$ci_level)) / 2)
diff_specs <- as.character(unlist(xr$difference_requires_specs))
scale <- cfg$axis$outcome_scale
out_name <- cfg$axis$outcome

ax <- lapply(contrast, function(re) readRDS(file.path(run_dir(re), "checkpoint",
                                                      "axis-designs.rds")))
names(ax) <- contrast

## One donor table over the union, with diagnosis, so the paired bootstrap can
## stratify exactly as the per-region bootstrap did. A donor's diagnosis is a
## donor attribute and must agree between regions.
dx_tab <- unique(rbindlist(lapply(contrast, function(re) {
    ck[[re]]$donors[, .(FID = as.character(FID), diagnosis)]
})))
if (anyDuplicated(dx_tab$FID)) stop("A donor carries different diagnoses in two regions")
dx_tab <- dx_tab[order(FID)]

region_estimate <- function(re, sp, sampled = NULL, drop_chrom = NULL) {
    des <- ck[[re]]$designs[[sp]]
    ids <- as.character(ck[[re]]$donors$FID[des$rows])
    Yr <- ck[[re]]$Y[des$rows, , drop = FALSE]
    idx <- if (is.null(sampled)) seq_along(ids) else {
        i <- match(sampled, ids); i[!is.na(i)]
    }
    X <- des$X[idx, , drop = FALSE]
    if (qr(X)$rank < ncol(X)) return(NA_real_)
    y <- age_outcomes(fit_age_matrix(X, Yr[idx, , drop = FALSE], des$age_col))[[out_name]]
    drop <- if (is.null(drop_chrom)) NULL else ck[[re]]$chrom == drop_chrom
    axis_estimate(ax[[re]][[sp]], y, scale, drop = drop)
}

draws <- list()
diff_rows <- list()
for (sp in diff_specs) {
    fitted_both <- all(vapply(contrast, function(re) !is.null(ck[[re]]$designs[[sp]]),
                              logical(1)))
    if (!fitted_both) stop("Spec ", sp, " was not fitted in both contrast regions")
    obs <- vapply(contrast, function(re) region_estimate(re, sp), numeric(1))
    d_obs <- obs[[1]] - obs[[2]]

    ## Donor half: paired bootstrap over the union, within diagnosis.
    bd <- matrix(NA_real_, nrow = B, ncol = 2, dimnames = list(NULL, contrast))
    for (b in seq_len(B)) {
        set.seed(seed_for(paste0("age-crossregion-", opts$cohort),
                          paste(contrast, collapse = "-"), sp, b))
        sampled <- dx_tab$FID[resample_rows(dx_tab$diagnosis)]
        for (re in contrast) bd[b, re] <- region_estimate(re, sp, sampled)
    }
    ok <- stats::complete.cases(bd)
    var_boot <- stats::var(bd[ok, 1] - bd[ok, 2])

    ## VMR half: JOINT delete-one-chromosome jackknife of the difference, the
    ## same chromosome removed from both regions at once, weighted by the
    ## combined VMR count per chromosome.
    chroms <- sort(unique(unlist(lapply(contrast, function(re) ck[[re]]$chrom))))
    n_tot <- sum(vapply(contrast, function(re) length(ck[[re]]$chrom), numeric(1)))
    hj <- n_tot / vapply(chroms, function(k) {
        sum(vapply(contrast, function(re) sum(ck[[re]]$chrom == k), numeric(1)))
    }, numeric(1))
    theta <- vapply(chroms, function(k) {
        region_estimate(contrast[1], sp, drop_chrom = k) -
            region_estimate(contrast[2], sp, drop_chrom = k)
    }, numeric(1))
    pseudo <- hj * d_obs - (hj - 1) * theta
    var_jk <- sum((pseudo - mean(pseudo))^2 / (hj - 1)) / length(chroms)

    se <- sqrt(var_boot + var_jk)
    ci <- d_obs + c(-1, 1) * z_ci * se
    diff_rows[[sp]] <- data.table(
        spec = sp, region_a = contrast[1], region_b = contrast[2],
        estimate_a = obs[[1]], estimate_b = obs[[2]],
        difference = d_obs, se = se,
        se_bootstrap = sqrt(var_boot), se_jackknife = sqrt(var_jk),
        ci_lower = ci[1], ci_upper = ci[2],
        ci_excludes_zero = ci[1] > 0 || ci[2] < 0,
        se_difference_if_independent = sqrt(stats::var(bd[ok, 1]) + stats::var(bd[ok, 2])),
        bootstrap_cor_between_regions = stats::cor(bd[ok, 1], bd[ok, 2]),
        n_bootstrap = B, n_bootstrap_failed = sum(!ok),
        n_chromosome_blocks = length(chroms),
        n_union_donors = nrow(dx_tab),
        n_shared_donors = length(Reduce(intersect, donors[contrast])),
        outcome = out_name, outcome_scale = scale,
        inference = "paired_donor_bootstrap_var_plus_joint_chromosome_jackknife_var")
    draws[[sp]] <- data.table(spec = sp, draw = seq_len(B),
                              estimate_a = bd[, 1], estimate_b = bd[, 2])
}
diffs <- rbindlist(diff_rows)
difference_claimed <- all(diffs$ci_excludes_zero)

## ------------------------------------------------------------- the decision
noncaudate <- per_region[!region %in% confounded]
token <- if (region_general) {
    if (difference_claimed) "REGION_GENERAL_WITH_DLPFC_HIPPOCAMPUS_DIFFERENCE" else
        "REGION_GENERAL"
} else if (sum(noncaudate$region_supported %in% TRUE) == 1L) {
    "SINGLE_NONCAUDATE_REGION"
} else if (any(noncaudate$region_reading == "OPPOSITE_DIRECTION")) {
    "OPPOSITE_DIRECTION"
} else {
    "NOT_SUPPORTED"
}
stale <- sum(!per_region$upstream_current)
smoke <- any(per_region$smoke_run %in% TRUE)

decision <- data.table(
    cohort = opts$cohort,
    aging_axis_association = token,
    n_regions_supported = n_supported,
    supported_regions = paste(per_region$region[supported], collapse = ","),
    supported_outside_caudate = sum(supported & outside),
    region_general = region_general,
    identified_difference_claimed = difference_claimed,
    identified_difference_requires_specs = paste(diff_specs, collapse = ","),
    caudate_tier = "descriptive_only",
    regions_are_independent_replicates = FALSE,
    max_pairwise_jaccard = max(overlap$jaccard),
    pooled_p_emitted = FALSE,
    manuscript_placement = as.character(cfg$interpretation$manuscript_placement),
    n_regions_citing_superseded_module_02 = stale,
    citable = !allow_unaccepted && stale == 0L && !smoke,
    built_with_unaccepted_runs = allow_unaccepted,
    runs = paste(paste0(regions, "=", runs), collapse = ","),
    cross_sectional_design = TRUE,
    causal_interpretation_allowed = FALSE,
    environmentally_determined_claim_allowed = FALSE,
    cross_region_raw_score_comparison_allowed = FALSE
)

sfx <- if (allow_unaccepted) paste0("-", opts$cohort, "-UNACCEPTED") else
    paste0("-", opts$cohort)
write_atomic(per_region, file.path(out_dir, paste0("aging-per-region", sfx, ".tsv")))
write_atomic(overlap, file.path(out_dir, paste0("aging-region-donor-overlap", sfx, ".tsv")))
write_atomic(diffs, file.path(out_dir, paste0("aging-identified-difference", sfx, ".tsv")))
write_atomic(rbindlist(draws), file.path(out_dir,
             paste0("aging-identified-difference-bootstrap", sfx, ".tsv")))
write_atomic(decision, file.path(out_dir, paste0("aging-cross-region-decision", sfx, ".tsv")))
if (nrow(ann_all)) {
    write_atomic(ann_all, file.path(out_dir, paste0("aging-annotation-associations", sfx, ".tsv")))
    write_atomic(ann_summary, file.path(out_dir,
                 paste0("aging-annotation-cross-region", sfx, ".tsv")))
}

print(per_region[, .(region, tier, region_reading, primary_estimate,
                     primary_p, upstream_current)])
print(overlap)
print(diffs[, .(spec, n_shared_donors, bootstrap_cor_between_regions,
                se, se_bootstrap, se_jackknife, se_difference_if_independent)])
message("[05] aging_axis_association = ", token,
        if (!decision$citable) "  (NOT CITABLE: unaccepted, smoke or stale inputs)" else "")
