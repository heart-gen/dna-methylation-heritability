#!/usr/bin/env Rscript
#### 09b_aging_application -- age-effect magnitude along the control axis ####
##
## Usage:
##   Rscript _h/02_axis_test.R --run-id age-AA-dlpfc-YYYYMMDD
##
## PRIMARY, per spec:
##   (beta_hat^2 - SE^2) / mean  ~  local_snp_contribution_score_z + technical covariates
##
## The outcome is the debiased squared age slope. Its expectation is beta^2
## whatever the SE, which is the whole point: the age model carries each VMR's
## local genetic variance in its residual, so |beta_hat|, its rank, |t| and
## -log10 p all couple to the score through the SE (config/aging.yml header
## records the design this replaced, and why). Dividing by the region mean
## makes the coefficient the proportional change in mean squared age effect per
## SD of score, which is what stage 05 compares between regions.
##
## Inference: SE^2 = donor-bootstrap variance + delete-one-chromosome block
## jackknife variance. The two halves are the donor-level and the VMR-level
## sampling; either alone under-covers (tests/). The bootstrap is used for its
## variance only, never a percentile interval: the bootstrap distribution of
## beta^2 - SE^2 is shifted by about SE^2 and would reintroduce the coupling.
##
## Arms change only the axis model (covariates or row subset) on the primary
## spec's age effects; the secondary outcome changes only the outcome.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
H_DIR <- Sys.getenv("V2_RUN_CODE",
                    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09b_aging_application", "_h"))
source(file.path(H_DIR, "run_config.R"))
source(file.path(H_DIR, "age_functions.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "09b_aging_application"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
if (run_is_sealed(manifest)) stop("Run is sealed and immutable: ", opts$run_id)
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mf("cohort"); region <- mf("region")
cfg <- load_run_config("aging", run_dir)
B <- as.integer(mf("n_bootstrap"))
if (!is.finite(B) || B < 2L) stop("Manifest carries no usable n_bootstrap")
predictor <- cfg$axis$predictor
base_covs <- as.character(unlist(cfg$axis$covariates))
alpha <- as.numeric(cfg$inference$alpha)

ck <- readRDS(file.path(run_dir, "checkpoint", "age-inputs.rds"))
if (!identical(ck$run_id, opts$run_id)) stop("Checkpoint belongs to another run")

## ------------------------------------------------------------- VMR features
score <- load_local_genetic_control(mf("upstream_local_genetic_variance_run_id"),
                                    region = region, cohort = cohort)
feat_f <- file.path(repo_root(), "04_repeat_repressive_architecture", "_m", "runs",
                    mf("upstream_repeat_architecture_run_id"), "results",
                    "vmr-features.tsv")
if (!file.exists(feat_f)) stop("Missing Module 04 features: ", feat_f)
feat <- fread(feat_f)
banned <- intersect(unlist(cfg$forbidden_columns), c(names(score), names(feat)))
if (length(banned)) {
    stop("Upstream table carries superseded column(s): ",
         paste(banned, collapse = ", "), " (AGENTS.md 3)")
}
if (!identical(unique(feat$vmr_set_id), mf("vmr_set_id"))) {
    stop("Module 04 feature table describes a different vmr_set_id")
}
arm_covs <- unique(unlist(lapply(cfg$axis$arms, function(a) a$add_covariates)))
arm_cols <- unique(unlist(lapply(cfg$axis$arms, function(a) a$subset$column)))
need <- unique(c(base_covs, arm_covs, arm_cols))
missing <- setdiff(need, names(feat))
if (length(missing)) {
    stop("Module 04 feature table lacks: ", paste(missing, collapse = ", "))
}
## The score comes from Module 02 directly, never from Module 04's copy, which
## may predate the accepted rescore. Only the named columns are taken, so the
## audit-only raw PVE estimate cannot reach a model.
dt <- data.table(vmr_id = ck$vmr_ids)
dt <- merge(dt, score[, c("vmr_id", predictor, "local_snp_contribution_quartile"),
                      with = FALSE], by = "vmr_id", all.x = TRUE, sort = FALSE)
dt <- merge(dt, feat[, c("vmr_id", need), with = FALSE], by = "vmr_id",
            all.x = TRUE, sort = FALSE)
dt <- dt[match(ck$vmr_ids, vmr_id)]
if (!identical(dt$vmr_id, ck$vmr_ids) || anyNA(dt[[predictor]])) {
    stop("Score/feature join lost or reordered VMRs")
}

## ------------------------------------------------------------- the test list
tests <- list()
add_test <- function(spec, arm, outcome, scale, role, covs, rows = NULL) {
    key <- paste(spec, arm, outcome, sep = "|")
    tests[[key]] <<- list(spec = spec, arm = arm, outcome = outcome, scale = scale,
                          role = role, ax = prepare_axis(dt, predictor, covs, rows))
}
spec_cfg <- cfg$age_model$specs
main_out <- cfg$axis$outcome; main_scale <- cfg$axis$outcome_scale
for (sp in names(ck$designs)) {
    add_test(sp, "base", main_out, main_scale, spec_cfg[[sp]]$role, base_covs)
}
for (arm in names(cfg$axis$arms)) {
    a <- cfg$axis$arms[[arm]]
    covs <- c(base_covs, as.character(unlist(a$add_covariates)))
    rows <- if (!is.null(a$subset)) {
        v <- dt[[a$subset$column]]
        is.finite(v) & v >= as.numeric(a$subset$min)
    } else NULL
    add_test("primary", arm, main_out, main_scale, a$role, covs, rows)
}
for (out in names(cfg$axis$secondary_outcomes)) {
    add_test("primary", "base", out, cfg$axis$secondary_outcomes[[out]]$scale,
             "secondary", base_covs)
}
by_spec <- split(names(tests), vapply(tests, `[[`, character(1), "spec"))

estimate_all <- function(outs, keys) {
    vapply(keys, function(k) {
        t <- tests[[k]]
        axis_estimate(t$ax, outs[[t$outcome]], t$scale)
    }, numeric(1))
}

## ------------------------------------------------------------- observed
Y <- ck$Y
obs_out <- lapply(ck$designs, function(des) {
    age_outcomes(fit_age_matrix(des$X, Y[des$rows, , drop = FALSE], des$age_col))
})
est <- unlist(lapply(names(by_spec), function(sp) estimate_all(obs_out[[sp]], by_spec[[sp]])))

## ------------------------------------------------------------- donor bootstrap
boot <- matrix(NA_real_, nrow = B, ncol = length(tests),
               dimnames = list(NULL, names(tests)))
message("[09b] ", B, " donor bootstrap draws x ", length(ck$designs), " specs x ",
        length(tests), " axis tests")
for (b in seq_len(B)) {
    for (sp in names(by_spec)) {
        des <- ck$designs[[sp]]
        set.seed(seed_for(opts$run_id, region, paste0("boot-", sp), b))
        i <- resample_rows(des$strata)
        X <- des$X[i, , drop = FALSE]
        if (qr(X)$rank < ncol(X)) next          # counted below, never imputed
        outs <- age_outcomes(fit_age_matrix(X, Y[des$rows, , drop = FALSE][i, , drop = FALSE],
                                            des$age_col))
        boot[b, by_spec[[sp]]] <- estimate_all(outs, by_spec[[sp]])
    }
    if (b %% 100L == 0L) message("[09b] bootstrap ", b, "/", B)
}

## ------------------------------------------------------------- inference
res <- rbindlist(lapply(names(tests), function(k) {
    t <- tests[[k]]
    y <- obs_out[[t$spec]][[t$outcome]]
    se_jk <- block_jackknife_se(t$ax, y, ck$chrom, t$scale)
    inf <- combined_inference(est[[k]], boot[, k], se_jk, alpha)
    data.table(
        test = k, spec = t$spec, arm = t$arm, outcome = t$outcome,
        outcome_scale = t$scale, role = t$role,
        n_vmrs = sum(t$ax$ok),
        n_chromosome_blocks = length(unique(ck$chrom[t$ax$ok])),
        estimate = est[[k]], se = inf$se,
        se_bootstrap = inf$se_bootstrap, se_jackknife = inf$se_jackknife,
        z = inf$z, p = inf$p, ci_lower = inf$ci_lower, ci_upper = inf$ci_upper,
        n_bootstrap = B, n_bootstrap_failed = B - inf$n_bootstrap_used,
        ## Diagnostic only: how far the bootstrap distribution sits from the
        ## estimate. Large for the debiased outcome by construction (~SE^2),
        ## which is why no percentile interval is formed.
        bootstrap_mean_minus_estimate = mean(boot[, k], na.rm = TRUE) - est[[k]],
        covariates = paste(t$ax$covariates, collapse = ","))
}))
res[, direction := fifelse(estimate < 0, "negative", "positive")]
res[, hypothesized_direction := fifelse(outcome == main_out,
                                        cfg$axis$hypothesized_sign, NA_character_)]
res[, `:=`(cohort = cohort, region = region, run_id = opts$run_id,
           predictor = predictor,
           inference = "donor_bootstrap_var_plus_chromosome_jackknife_var")]

## ------------------------------------------------------------- descriptive
d <- obs_out$primary[[main_out]]
q <- data.table(quartile = dt$local_snp_contribution_quartile, d = d,
                signed = obs_out$primary$signed_beta_age)
quart <- q[, .(n_vmrs = .N,
               mean_debiased_sq_beta_relative = mean(d) / mean(q$d),
               rms_age_effect_per_decade = sqrt(max(mean(d), 0)),
               frac_gain_with_age = mean(signed > 0)),
           by = quartile][order(quartile)]
quart[, `:=`(region = region, run_id = opts$run_id, spec = "primary",
             descriptive_only = TRUE)]

write_atomic(res, file.path(run_dir, "results", "axis-tests.tsv"))
write_atomic(quart, file.path(run_dir, "results", "axis-quartile-summary.tsv"))
saveRDS(boot, file.path(run_dir, "checkpoint", "donor-bootstrap.rds"))
## Each spec's base axis design, for stage 05's paired bootstrap.
base_keys <- names(tests)[vapply(tests, function(t) {
    t$arm == "base" && t$outcome == main_out
}, logical(1))]
saveRDS(stats::setNames(lapply(base_keys, function(k) tests[[k]]$ax),
                        vapply(tests[base_keys], `[[`, character(1), "spec")),
        file.path(run_dir, "checkpoint", "axis-designs.rds"))
append_manifest(list(dir = run_dir), list(
    n_axis_tests = as.character(nrow(res)),
    n_bootstrap_completed = as.character(B),
    max_bootstrap_failed = as.character(max(res$n_bootstrap_failed))
))
print(res[, .(spec, arm, outcome, n_vmrs, estimate = signif(estimate, 3),
              se = signif(se, 3), se_boot = signif(se_bootstrap, 3),
              se_jk = signif(se_jackknife, 3), p = signif(p, 3))])
