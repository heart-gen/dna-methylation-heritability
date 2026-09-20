#### 10_environmental_exploratory -- one definition of the per-VMR model ####
##
## Stage 2 fits the per-VMR exposure model; stage 3's donor bootstrap refits the
## SAME model thousands of times on resampled donors. Two copies of the stratum
## filter, the covariate list or the design matrix would be a silent divergence
## waiting to happen -- the bootstrap would then be resampling a model that is
## not the one whose estimate it is putting a SE on. Both stages call these.
##
## Added 2026-09-19 with the debiased outcome and the bootstrap (PI approval of
## the same date; see config/environmental.yml:testing).

## Covariates for one stratum. Inside a single-diagnosis stratum primarydx is
## constant, so keeping it would make the design rank-deficient; dropping it is
## not a weaker adjustment, because the collider path it was adjusting for is
## closed by the restriction itself.
exposure_base_terms <- function(st, env, pc_names = NULL) {
    keep_dx <- isTRUE(env$strata[[st]]$include_diagnosis_covariate)
    c("age", "factor(sex)",
      if (keep_dx) "factor(diagnosis)" else NULL,
      pc_names)
}

## Donors contributing to one exposure x stratum family. `md` is the merged
## donor table (phenotype + covariates + exposures) for one VMR.
## `phenotype` is separate because stage 3 works from a donor x VMR matrix and
## has no single phenotype column; it passes NULL and screens Y's finiteness
## per VMR itself. Stage 2 uses the default and so is unchanged.
exposure_keep_rows <- function(md, ex, st, env, phenotype = md$phenotype) {
    in_stratum <- if (identical(st, "all")) rep(TRUE, nrow(md)) else
        md$primarydx %in% as.character(unlist(env$strata[[st]]$dx_filter))
    keep <- in_stratum & !is.na(md[[ex]]) &
        is.finite(as.numeric(md$age)) & !is.na(md$sex) & !is.na(md$diagnosis)
    if (!is.null(phenotype)) keep <- keep & is.finite(as.numeric(phenotype))
    keep
}

exposure_formula <- function(ex, st, env, pc_names = NULL, response = ".y") {
    stats::reformulate(c(ex, exposure_base_terms(st, env, pc_names)),
                       response = response)
}

## The model matrix for one family, plus the column indices the exposure
## occupies (1 column for a binary exposure, k-1 for a k-level categorical one).
## `d` must already be restricted to the family's donors.
exposure_design <- function(d, ex, st, env, pc_names = NULL) {
    if (is.factor(d[[ex]])) d[[ex]] <- droplevels(d[[ex]])
    f <- exposure_formula(ex, st, env, pc_names, response = NULL)
    X <- stats::model.matrix(f, data = d)
    cols <- which(startsWith(colnames(X), ex))
    if (!length(cols)) stop("Exposure '", ex, "' contributes no column")
    list(X = X, cols = cols, df_term = length(cols))
}
