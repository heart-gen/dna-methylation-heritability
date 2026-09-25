#!/usr/bin/env Rscript
#### 07_transcription_splicing_coupling -- the coupling tests ####
##
## Usage:
##   Rscript _h/03_test_coupling.R --run-id tsc-AA-caudate-YYYYMMDD
##
## Three prespecified tests per modality, asking whether VMRs that are
## genetically regulated are more often transcriptionally coupled:
##
##   (a) any CpG meQTL support        -> coupling
##   (b) proportion of meQTL-supported CpGs -> coupling
##   (c) continuous local SNP contribution score -> coupling
##
## Ported from the legacy 01_meqtl_tx_enrichment.py with the predictors
## swapped. The legacy version regressed on `r_squared_cv`, `h2_category` and
## `h2_unscaled`; all three are banned (AGENTS.md 3), and the guard below stops
## the run if any of them reaches this stage.
##
## Allowed reading of a positive result: "genetically regulated VMRs are more
## frequently transcriptionally coupled." NOT that methylation mediates a
## genetic effect on expression or splicing -- nothing here tests mediation.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE", file.path(Sys.getenv("V2_REPO_ROOT", "."), "07_transcription_splicing_coupling", "_h")), "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(sandwich)
    library(lmtest)
})

MODULE <- "07_transcription_splicing_coupling"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field, required = TRUE) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0 || is.na(v[1])) {
        if (required) stop("Manifest has no value for '", field, "'")
        return(NA_character_)
    }
    v[1]
}
cohort <- mf("cohort"); region <- mf("region")
burden_run <- mf("upstream_cpg_meqtl_burden_run_id")
ts <- load_run_config("transcription_splicing", run_dir)
enabled <- strsplit(mf("modalities"), ",", fixed = TRUE)[[1]]

## ------------------------------------------------------------ predictors
burden_f <- file.path(repo_root(), "05_cpg_meqtl_burden", "_m", "runs",
                      burden_run, "results", "vmr-meqtl-burden.tsv")
if (!file.exists(burden_f)) stop("Module 05 burden table not found: ", burden_f)
burden <- fread(burden_f)

## Fail closed on any superseded estimator column. The legacy link tables carry
## all of these; a table that still has one is a pre-repair artifact and must
## not reach a model (AGENTS.md 3).
banned <- intersect(ts$forbidden_columns, names(burden))
if (length(banned)) {
    stop("Module 05 burden table carries banned column(s): ",
         paste(banned, collapse = ", "))
}
needed <- c("vmr_id", "local_snp_contribution_score_z",
            "proportion_cpgs_with_sig_meqtl", "n_cpgs_with_sig_meqtl",
            "n_tested_cpgs", "vmr_length", "cpg_density")
missing <- setdiff(needed, names(burden))
if (length(missing)) {
    stop("Module 05 burden table is missing: ", paste(missing, collapse = ", "))
}

burden[, any_meqtl_support := as.integer(n_cpgs_with_sig_meqtl > 0)]

## num_snps is a required adjustment (AGENTS.md 7.6) and lives in the Module 02
## score table rather than the burden table.
lgv_run <- mf("upstream_local_genetic_variance_run_id")
lgv <- load_local_genetic_control(lgv_run, region = region, cohort = cohort,
                                  eligible_only = FALSE)
burden <- merge(burden,
                unique(lgv[, .(vmr_id, num_snps, methylation_variance)],
                       by = "vmr_id"),
                by = "vmr_id", all.x = TRUE)

results <- list()
n_vmrs_modality <- list()
per_modality_tables <- list()

for (mod in enabled) {
    sum_f <- file.path(run_dir, "results", paste0(mod, "-vmr-summary.tsv"))
    if (!file.exists(sum_f)) {
        warning("No VMR summary for ", mod, "; skipping")
        next
    }
    vs <- fread(sum_f)
    d <- merge(vs, burden, by = "vmr_id")
    ## A zero-row join means the two tables key their VMRs differently, not that
    ## there is nothing to test. That happened once already (module 07 emitted
    ## the legacy underscore form while 02/05 use chr:start-end) and surfaced
    ## only as "insufficient variation" with n = 0 in every test, which reads
    ## like a null result. Fail loudly instead.
    if (nrow(d) == 0) {
        stop(mod, ": joining the VMR summary to the Module 05 burden table on ",
             "vmr_id produced zero rows. Example summary id '", vs$vmr_id[1],
             "', example burden id '", burden$vmr_id[1],
             "'. These are different identifier conventions, not an empty result.")
    }
    ## The locked power floor. `config/transcription_splicing.yml:gates` declares
    ## min_vmrs_tested, and a modality below it is retained and reported but kept
    ## OUT of the coupling-test FDR family -- see the note at the p.adjust call.
    ## This used to warn and nothing more, while 04_apply_gates.R evaluated the
    ## same floor against max() across modalities, so the strongest modality
    ## satisfied it for every other one and the floor never bound anything.
    n_vmrs_modality[[mod]] <- nrow(d)
    if (nrow(d) < ts$gates$min_vmrs_tested) {
        warning(mod, ": only ", nrow(d), " VMRs after joining predictors, ",
                "below the locked min_vmrs_tested (", ts$gates$min_vmrs_tested,
                "); its tests are reported but excluded from the FDR family")
    }

    banned2 <- intersect(ts$forbidden_columns, names(d))
    if (length(banned2)) {
        stop(mod, ": banned column(s) reached the model frame: ",
             paste(banned2, collapse = ", "))
    }

    d[, coupled := as.integer(any_sig_fdr > 0)]
    ## log1p on the count-like adjustments: they are heavily right-skewed, and
    ## an untransformed feature count would let a handful of gene-dense loci
    ## dominate the fit.
    d[, `:=`(log_n_features = log1p(n_features_tested),
             log_vmr_length = log1p(vmr_length),
             log_num_snps = log1p(num_snps),
             log_min_distance = log1p(min_distance_to_feature),
             log_n_tested_cpgs = log1p(n_tested_cpgs))]

    covar_terms <- c("log_n_features", "log_vmr_length", "log_min_distance",
                     "methylation_variance", "log_num_snps", "cpg_density",
                     "log_n_tested_cpgs")
    covar_terms <- covar_terms[vapply(covar_terms, function(cc) {
        cc %in% names(d) && sum(is.finite(d[[cc]])) > 0 &&
            data.table::uniqueN(d[[cc]]) > 1
    }, logical(1))]

    for (pname in names(ts$coupling$predictors)) {
        pcol <- switch(pname,
                       any_meqtl_support = "any_meqtl_support",
                       meqtl_proportion = "proportion_cpgs_with_sig_meqtl",
                       local_genetic_control = "local_snp_contribution_score_z")
        if (!pcol %in% names(d)) {
            warning(mod, "/", pname, ": predictor ", pcol, " absent")
            next
        }
        keep <- complete.cases(d[, c("coupled", pcol, covar_terms), with = FALSE])
        dd <- d[keep]
        if (nrow(dd) < 50 || uniqueN(dd$coupled) < 2 ||
            uniqueN(dd[[pcol]]) < 2) {
            results[[length(results) + 1]] <- data.table(
                modality = mod, predictor = pname, predictor_column = pcol,
                n = nrow(dd), estimate = NA_real_, se = NA_real_,
                z = NA_real_, p = NA_real_,
                note = "insufficient variation")
            next
        }
        form <- stats::as.formula(paste("coupled ~", pcol, "+",
                                        paste(covar_terms, collapse = " + ")))
        ## Binary outcome over loci with substantial overdispersion across the
        ## genome; quasibinomial with HC3 sandwich SEs matches the Module 04
        ## convention for the same kind of per-locus outcome.
        fit <- try(stats::glm(form, data = dd,
                              family = stats::quasibinomial(link = "logit")),
                   silent = TRUE)
        if (inherits(fit, "try-error")) {
            results[[length(results) + 1]] <- data.table(
                modality = mod, predictor = pname, predictor_column = pcol,
                n = nrow(dd), estimate = NA_real_, se = NA_real_,
                z = NA_real_, p = NA_real_, note = "model failed")
            next
        }
        ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
        if (!pcol %in% rownames(ct)) next
        results[[length(results) + 1]] <- data.table(
            modality = mod, predictor = pname, predictor_column = pcol,
            n = nrow(dd),
            n_coupled = sum(dd$coupled),
            estimate = ct[pcol, 1], se = ct[pcol, 2],
            z = ct[pcol, 3], p = ct[pcol, 4],
            covariates = paste(covar_terms, collapse = "+"),
            note = "")
    }
    per_modality_tables[[mod]] <- d[, .(vmr_id, modality, coupled,
                                        n_features_tested, n_sig_fdr,
                                        any_meqtl_support,
                                        proportion_cpgs_with_sig_meqtl,
                                        local_snp_contribution_score_z)]
}

if (length(results) == 0) stop("No coupling test could be fitted")
res <- rbindlist(results, fill = TRUE)
res[, `:=`(cohort = cohort, region = region, run_id = opts$run_id)]

## ------------------------------------------------- the coupling-test FDR family
##
## Correct across the tests actually reported -- predictors x modalities -- but
## only for modalities that clear the locked power floor.
##
## `expression_abc` links 243-305 VMRs against a declared min_vmrs_tested of 500,
## and with roughly 8-15 of those coupled a binary predictor has no off-cell to
## estimate from: `any_meqtl_support` separated completely in all three AA regions
## on 2026-09-25, returning a boundary estimate near 17.6 with an SE near 0.4 and
## p underflowing to 0. BH then sorted that unidentified coefficient first, so the
## most significant-looking row in the module was a fit that had not converged.
##
## The exclusion is NOT a reaction to that result. It applies the floor the PI
## locked in config, evaluated per modality as the key is written, and it is
## decided on the modality's SIZE -- which is a property of the ABC link set, not
## of any p-value. An underpowered modality is retained in the table with a
## nominal p and no q, the pattern 06_partitioned_heritability uses for its
## membership tau, so the reader sees the estimate and the reason rather than a
## silent gap.
.nv <- unlist(n_vmrs_modality)
res[, n_vmrs_modality := as.integer(.nv[modality])]
## A modality with no recorded size must not become NA and slip through the
## `any()` below as neither in nor out of the family; absent means excluded.
if (anyNA(res$n_vmrs_modality)) {
    warning("no recorded VMR count for modality(s) ",
            paste(unique(res$modality[is.na(res$n_vmrs_modality)]),
                  collapse = ", "), "; treating as below the power floor")
}
res[, power_eligible := !is.na(n_vmrs_modality) &
        n_vmrs_modality >= ts$gates$min_vmrs_tested]
res[, in_fdr_family := power_eligible & is.finite(p)]
res[, fdr_exclusion_reason := fifelse(
        power_eligible, NA_character_,
        sprintf("modality links %d VMRs, below the locked min_vmrs_tested of %d",
                n_vmrs_modality, ts$gates$min_vmrs_tested))]
res[, q := NA_real_]
if (any(res$in_fdr_family)) {
    res[in_fdr_family == TRUE,
        q := p.adjust(p, method = ts$association$fdr_method)]
}
res[, fdr_family_size := sum(res$in_fdr_family)]
setorder(res, p)
write_atomic(res, file.path(run_dir, "results", "coupling-tests.tsv"))

if (length(per_modality_tables)) {
    write_atomic(rbindlist(per_modality_tables, fill = TRUE),
                 file.path(run_dir, "results", "coupling-model-frame.tsv"))
}

print(res[, .(modality, predictor, n, n_coupled, estimate, p, q)])
message("[07] coupling tests written")
