#!/usr/bin/env Rscript
#### 10_environmental_exploratory -- per-VMR exposure association (one chrom) ####
##
## Usage:
##   Rscript _h/02_vmr_exposure_association.R --run-id env-... --chrom 7
##
## Stage A of the module, and the direct v2 descendant of v1's
## environmental-analysis/*/correlation/_h/01.corr_pheno.R. Two things changed.
##
## The locus universe is Module 02's ELIGIBLE loci, so the exposure scan and the
## local-control score describe the same set and the Stage B join has no
## silently missing rows.
##
## The covariate model comes from 00_shared/locus_io.R rather than being written
## out here. v1 used `age + sex + dx + afr_ances`; the v2 arms residualize
## methylation on pooled snpPC1-3 during VMR discovery, so a global ancestry
## proportion is both superseded and not the same adjustment. Reusing the shared
## reader means an exposure association conditions on exactly what a local
## genetic-control estimate conditions on.
##
## No FDR here. This runs per chromosome; correcting inside a task would give 22
## families per exposure (AGENTS.md 10.3).

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "10_environmental_exploratory", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "10_environmental_exploratory"
opts <- parse_v2_args(require = c("run_id", "chrom"))
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
chrom_label <- sub("^chr", "", as.character(opts$chrom), ignore.case = TRUE)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mf("cohort"); region <- mf("region")
catalog_cohort <- if (is.na(mf("catalog_cohort"))) cohort else mf("catalog_cohort")
estimation_group <- if (is.na(mf("estimation_group"))) cohort else mf("estimation_group")
env <- load_run_config("environmental", run_dir)

## The unit is the exposure x stratum pair, encoded "exposure@stratum" in the
## manifest. Each pair is fitted separately and becomes its own BH family.
pairs_raw <- strsplit(mf("eligible_exposures"), ",", fixed = TRUE)[[1]]
pairs_raw <- pairs_raw[nzchar(pairs_raw)]
if (!length(pairs_raw)) stop("Manifest lists no eligible exposures; run stage 01.")
pairs <- data.table(
    exposure = sub("@.*$", "", pairs_raw),
    stratum = sub("^.*@", "", pairs_raw))
eligible <- unique(pairs$exposure)

expo <- fread(file.path(run_dir, "results", "exposure-matrix.tsv"),
              na.strings = c("NA", ""))
for (nm in names(env$exposures$categorical)) {
    if (nm %in% names(expo)) {
        expo[[nm]] <- factor(expo[[nm]],
                             levels = names(env$exposures$categorical[[nm]]$collapse))
    }
}

## Locus universe: Module 02's eligible loci on this chromosome.
lgv_run <- mf("upstream_local_genetic_variance_run_id")
if (is.na(lgv_run)) stop("Manifest carries no upstream_local_genetic_variance_run_id")
score <- load_local_genetic_control(lgv_run, region = region, cohort = cohort)
banned <- intersect(unlist(env$forbidden_columns), names(score))
if (length(banned)) {
    stop("Module 02 table carries superseded estimator column(s): ",
         paste(banned, collapse = ", "), " (AGENTS.md 3)")
}
audit <- intersect(unlist(env$audit_only_columns), names(score))
if (length(audit)) {
    ## AGENTS.md 7.2 keeps the raw estimate for audit inside Module 02. Dropping
    ## it here means no stage of this module can accidentally model or report it.
    score[, (audit) := NULL]
}
## `..chrom_label` is column-selection syntax and does not resolve in `i`, so
## the filter value carries a name that differs from the column's.
chrom_filter <- chrom_label
score[, chrom_label := sub("^chr", "", as.character(chrom), ignore.case = TRUE)]
tasks <- score[chrom_label == chrom_filter,
               .(vmr_id, chrom, start, end, n_cpgs)]
if (!nrow(tasks)) {
    message("[10] chr", chrom_label, ": no eligible VMRs")
    write_atomic(data.table(), file.path(run_dir, "results", "per-chrom",
                                         paste0("assoc-chr", chrom_label, ".tsv")))
    quit(save = "no", status = 0)
}

vmr_run_dir <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs",
                         mf("upstream_vmr_catalog_run_id"))
covar_prefix <- if (identical(catalog_cohort, "AA")) "TOPMed_LIBD.AA" else "TOPMed_LIBD"

fit_one <- function(task) {
    ph <- try(load_locus_phenotype(task = task, vmr_run_dir = vmr_run_dir,
                                   cohort = catalog_cohort,
                                   covar_prefix = covar_prefix), silent = TRUE)
    if (inherits(ph, "try-error")) {
        return(data.table(vmr_id = task$vmr_id, exposure = NA_character_,
                          stratum = NA_character_,
                          term = NA_character_, status = "phenotype_unavailable",
                          message = as.character(ph)))
    }
    md <- Reduce(function(x, y) merge(x, y, by = c("FID", "IID"), all = FALSE),
                 list(ph$phenotype, ph$covar, ph$qcovar))
    if (!is.null(ph$pcs)) {
        md <- merge(md, ph$pcs, by = c("FID", "IID"), all = FALSE)
    }
    md <- merge(md, expo, by.x = "FID", by.y = "brnum", all = FALSE)

    rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
        ex <- pairs$exposure[i]; st <- pairs$stratum[i]
        ## Inside a single-diagnosis stratum primarydx is constant, so keeping it
        ## as a covariate would make the design matrix rank-deficient. Dropping
        ## it is not a weaker adjustment: the collider path it was adjusting for
        ## is closed by the restriction itself.
        keep_dx <- isTRUE(env$strata[[st]]$include_diagnosis_covariate)
        base_terms <- c("age", "factor(sex)",
                        if (keep_dx) "factor(diagnosis)" else NULL,
                        ph$pc_names)
        in_stratum <- if (identical(st, "all")) rep(TRUE, nrow(md)) else
            md$primarydx %in% as.character(unlist(env$strata[[st]]$dx_filter))
        keep <- in_stratum & !is.na(md[[ex]]) &
            is.finite(as.numeric(md$phenotype)) &
            is.finite(as.numeric(md$age)) & !is.na(md$sex) & !is.na(md$diagnosis)
        ## md is a data.frame (the shared reader returns data.frames), so a
        ## single-argument subscript would select COLUMNS. Be explicit.
        d <- md[keep, , drop = FALSE]
        if (is.factor(d[[ex]])) d[[ex]] <- droplevels(d[[ex]])
        n_used <- nrow(d)
        classes <- table(d[[ex]])
        if (n_used < 3L || length(classes) < 2L) {
            return(data.table(vmr_id = task$vmr_id, exposure = ex, stratum = st,
                              term = NA_character_, status = "insufficient_variation",
                              n_used = n_used))
        }
        d$.y <- as.numeric(d$phenotype)
        f <- stats::reformulate(c(ex, base_terms), response = ".y")
        fit <- try(stats::lm(f, data = d), silent = TRUE)
        if (inherits(fit, "try-error")) {
            return(data.table(vmr_id = task$vmr_id, exposure = ex, stratum = st,
                              term = NA_character_, status = "model_failed",
                              n_used = n_used, message = as.character(fit)))
        }
        co <- summary(fit)$coefficients
        rn <- rownames(co)
        keep_terms <- startsWith(rn, ex)
        if (!any(keep_terms)) {
            return(data.table(vmr_id = task$vmr_id, exposure = ex, stratum = st,
                              term = NA_character_, status = "term_dropped",
                              n_used = n_used))
        }
        ## An exposure with k levels contributes k-1 terms. The per-VMR p-value
        ## carried into Stage B is the joint F test against the covariate-only
        ## model, not the smallest coefficient p -- taking a minimum across
        ## levels would be an undeclared selection step.
        null_fit <- stats::lm(stats::reformulate(base_terms, response = ".y"),
                              data = d)
        av <- stats::anova(null_fit, fit)
        joint_p <- av[["Pr(>F)"]][2]
        data.table(
            vmr_id = task$vmr_id, exposure = ex, stratum = st,
            term = rn[keep_terms],
            beta = co[keep_terms, "Estimate"],
            se = co[keep_terms, "Std. Error"],
            t = co[keep_terms, "t value"],
            p_term = co[keep_terms, "Pr(>|t|)"],
            p_joint = joint_p,
            df_num = av[["Df"]][2],
            n_used = n_used,
            status = "ok"
        )
    }), fill = TRUE)
}

res <- rbindlist(lapply(seq_len(nrow(tasks)), function(i) {
    fit_one(as.list(tasks[i]))
}), fill = TRUE)

res <- merge(res, tasks, by = "vmr_id", all.x = TRUE, sort = FALSE)
res[, `:=`(cohort = cohort, region = region, run_id = opts$run_id,
           estimation_group = estimation_group)]
write_atomic(res, file.path(run_dir, "results", "per-chrom",
                            paste0("assoc-chr", chrom_label, ".tsv")))
message("[10] chr", chrom_label, ": ", nrow(tasks), " VMRs x ",
        nrow(pairs), " exposure-stratum pairs -> ", nrow(res), " rows")
