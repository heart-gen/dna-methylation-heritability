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
##
## A GATING arm built from the DNAm scMD proportions is not fitted where scMD
## fails its integration gate, exactly as the `cell_scmd` age spec is not fitted
## there (01_age_effects.R). `cell_composition_r2` is such an arm: Module 04
## builds it from dnam-scmd-proportions-{region}.tsv, so gating a region's
## verdict on it while declining to fit cell_scmd in that same region would
## admit the same quantity under a second name. See age_functions.R:
## SCMD_DERIVED_FEATURES for the judgement and what would change if it is
## overturned. Declined arms go to results/axis-arms-skipped.tsv.
##
## DESCRIPTIVE annotation associations (config/aging.yml:
## annotation_associations) are fitted here too, because they need the same
## donor-bootstrap draws: the same outcomes regressed on one Module 04 or 07
## annotation at a time instead of the score. They are written to their own
## table and never reach the region reading.

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
## The scMD integration gate, re-read from the concordance table exactly as
## stage 01 read it, and cross-checked against what stage 01 recorded: the arm
## rule below and the `cell_scmd` spec rule must never disagree within a run.
scmd_ok <- scmd_gate_passes(region, cfg$cell_composition$scmd_gate_table)
gate_recorded <- mf("scmd_integration_gate")
if (!is.na(gate_recorded) &&
    !identical(gate_recorded, if (scmd_ok) "PASS" else "FAIL")) {
    stop("scMD integration gate is ", if (scmd_ok) "PASS" else "FAIL",
         " here but stage 01 recorded ", gate_recorded)
}

arm_scmd <- list()
arm_skipped <- list()
for (arm in names(cfg$axis$arms)) {
    a <- cfg$axis$arms[[arm]]
    arm_scmd[[arm]] <- arm_is_scmd_derived(a)
    ## A gating arm built from the scMD proportions is not fitted where scMD
    ## fails its integration gate -- the same treatment the `cell_scmd` spec
    ## gets in 01_age_effects.R, for the same reason (age_functions.R:
    ## SCMD_DERIVED_FEATURES documents why a derived scalar is still an scMD
    ## adjustment when it is used to gate).
    if (arm_skipped_for_scmd(a, scmd_ok)) {
        arm_skipped[[arm]] <- data.table(
            arm = arm, role = as.character(a$role), fitted = FALSE,
            reason = "scmd_integration_gate_fails_in_region",
            scmd_derived_covariate = TRUE,
            covariates = paste(as.character(unlist(a$add_covariates)),
                               collapse = ","))
        message("[09b] arm ", arm, " not fitted: scMD integration gate FAILS in ",
                region)
        next
    }
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

## ------------------------------------------------------------- annotation tests
## One indicator (or z-scored continuous annotation) at a time in place of the
## score, on the primary spec's age effects. Built here so the bootstrap loop
## below evaluates them on the very draws the axis tests use.
ann_cfg <- cfg$annotation_associations
ann_tests <- list()
ann_skipped <- list()
if (!is.null(ann_cfg)) {
    m04 <- ann_cfg$module_04
    bin_cols <- as.character(unlist(c(m04$chromatin_any, m04$repeat_any)))
    cont_cols <- as.character(unlist(m04$continuous_z))
    miss <- setdiff(c(bin_cols, cont_cols, "broad_genomic_annotation"), names(feat))
    if (length(miss)) stop("Module 04 feature table lacks: ", paste(miss, collapse = ", "))
    src <- feat[match(dt$vmr_id, feat$vmr_id)]
    ann_dt <- copy(dt)
    ann_src <- list()      # annotation column -> source module
    ann_extra <- list()    # annotation column -> extra covariates
    for (a in bin_cols) {
        ann_dt[[a]] <- as.numeric(as.logical(src[[a]]))
        ann_src[[a]] <- "04_repeat_repressive_architecture"
    }
    for (lv in as.character(unlist(m04$genomic_context_levels))) {
        if (!lv %in% src$broad_genomic_annotation) {
            stop("broad_genomic_annotation has no level '", lv, "'")
        }
        a <- paste0("genomic_", lv)
        ann_dt[[a]] <- as.numeric(src$broad_genomic_annotation == lv)
        ann_src[[a]] <- "04_repeat_repressive_architecture"
    }
    for (cc in cont_cols) {
        a <- paste0(cc, "_z")
        ann_dt[[a]] <- as.numeric(scale(src[[cc]]))
        ann_src[[a]] <- "04_repeat_repressive_architecture"
    }
    ## Module 07: the per-VMR coupled indicator, adjusted for how many features
    ## the VMR was tested against.
    m07 <- ann_cfg$module_07
    tsc_f <- file.path(repo_root(), "07_transcription_splicing_coupling", "_m", "runs",
                       mf("upstream_transcription_coupling_run_id"), "results",
                       "coupling-model-frame.tsv")
    if (!file.exists(tsc_f)) stop("Missing Module 07 model frame: ", tsc_f)
    tsc <- fread(tsc_f)
    for (mod in as.character(unlist(m07$modalities))) {
        m <- tsc[modality == mod]
        a <- paste0("coupled_", mod); nf <- paste0("n_features_tested_", mod)
        idx <- match(dt$vmr_id, m$vmr_id)
        ann_dt[[a]] <- as.numeric(m$coupled[idx])
        ann_dt[[nf]] <- as.numeric(m$n_features_tested[idx])
        ann_src[[a]] <- "07_transcription_splicing_coupling"
        ann_extra[[a]] <- nf
    }
    min_ann <- as.integer(ann_cfg$min_annotated_vmrs)
    min_cpl <- as.integer(m07$min_coupled_vmrs)
    for (a in names(ann_src)) {
        x <- ann_dt[[a]]
        binary <- all(x[is.finite(x)] %in% c(0, 1))
        n_ann <- if (binary) sum(x == 1, na.rm = TRUE) else NA_integer_
        n_oth <- if (binary) sum(x == 0, na.rm = TRUE) else NA_integer_
        floor_n <- if (identical(ann_src[[a]], "07_transcription_splicing_coupling")) min_cpl else min_ann
        if (binary && min(n_ann, n_oth) < floor_n) {
            ann_skipped[[a]] <- data.table(annotation = a, source_module = ann_src[[a]],
                n_annotated = n_ann, reason = paste0("fewer than ", floor_n,
                                                     " VMRs in a class"))
            next
        }
        for (adj in as.character(unlist(ann_cfg$adjustments))) {
            covs <- c(base_covs, if (adj == "technical_score") predictor,
                      unlist(ann_extra[[a]]))
            ax <- prepare_axis(ann_dt, a, covs)
            for (out in names(ann_cfg$outcomes)) {
                ann_tests[[paste(a, adj, out, sep = "|")]] <- list(
                    annotation = a, source_module = ann_src[[a]], adjustment = adj,
                    outcome = out, scale = as.character(ann_cfg$outcomes[[out]]),
                    binary = binary, ax = ax)
            }
        }
    }
    if (!identical(as.character(ann_cfg$spec), "primary")) {
        stop("annotation_associations are prespecified on the primary spec only")
    }
}
estimate_ann <- function(outs) {
    vapply(ann_tests, function(t) axis_estimate(t$ax, outs[[t$outcome]], t$scale),
           numeric(1))
}

scmd_arms <- names(arm_scmd)[vapply(arm_scmd, isTRUE, logical(1))]

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
est_ann <- estimate_ann(obs_out$primary)

## ------------------------------------------------------------- donor bootstrap
boot <- matrix(NA_real_, nrow = B, ncol = length(tests),
               dimnames = list(NULL, names(tests)))
boot_ann <- matrix(NA_real_, nrow = B, ncol = length(ann_tests),
                   dimnames = list(NULL, names(ann_tests)))
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
        if (sp == "primary" && length(ann_tests)) boot_ann[b, ] <- estimate_ann(outs)
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
        ## Whether this row's adjustment comes from the DNAm scMD proportions,
        ## so a reader can see which rows the integration gate governs without
        ## knowing how cell_composition_r2 is built.
        scmd_derived_covariate = t$arm %in% scmd_arms ||
            isTRUE(spec_cfg[[t$spec]]$requires_scmd_integration_gate),
        covariates = paste(t$ax$covariates, collapse = ","))
}))
res[, direction := fifelse(estimate < 0, "negative", "positive")]
res[, hypothesized_direction := fifelse(outcome == main_out,
                                        cfg$axis$hypothesized_sign, NA_character_)]
res[, `:=`(cohort = cohort, region = region, run_id = opts$run_id,
           predictor = predictor,
           scmd_integration_gate = if (scmd_ok) "PASS" else "FAIL",
           inference = "donor_bootstrap_var_plus_chromosome_jackknife_var")]

## Arms declined in this region, written even when empty so stage 03 always has
## a reason to record for a member it cannot find.
skipped <- if (length(arm_skipped)) rbindlist(arm_skipped, fill = TRUE) else
    data.table(arm = character(0), role = character(0), fitted = logical(0),
               reason = character(0), scmd_derived_covariate = logical(0),
               covariates = character(0))
skipped[, `:=`(cohort = cohort, region = region, run_id = opts$run_id,
               scmd_integration_gate = if (scmd_ok) "PASS" else "FAIL")]

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

## ------------------------------------------------------------- annotation table
ann_res <- rbindlist(lapply(names(ann_tests), function(k) {
    t <- ann_tests[[k]]
    y <- obs_out$primary[[t$outcome]]
    inf <- combined_inference(est_ann[[k]], boot_ann[, k],
                              block_jackknife_se(t$ax, y, ck$chrom, t$scale), alpha)
    x <- ann_dt[[t$annotation]][t$ax$ok]
    yy <- y[t$ax$ok]
    if (t$scale == "relative_to_mean") yy <- yy / mean(yy)
    g <- obs_out$primary$signed_beta_age[t$ax$ok]
    data.table(
        annotation = t$annotation, source_module = t$source_module,
        adjustment = t$adjustment, outcome = t$outcome, outcome_scale = t$scale,
        annotation_type = if (t$binary) "indicator" else "continuous_z",
        ## Descriptive rows are fitted in every region on purpose
        ## (config/aging.yml:213 keeps cell_composition_r2 here with the caveat
        ## that outside caudate it is a weak composition proxy). Flagged so a
        ## reader of this table cannot mistake it for a gated composition
        ## adjustment; it never enters the region reading.
        scmd_derived_annotation = sub("_z$", "", t$annotation) %in% SCMD_DERIVED_FEATURES,
        n_vmrs = sum(t$ax$ok),
        n_annotated = if (t$binary) sum(x == 1) else NA_integer_,
        estimate = est_ann[[k]], se = inf$se,
        se_bootstrap = inf$se_bootstrap, se_jackknife = inf$se_jackknife,
        z = inf$z, p = inf$p, ci_lower = inf$ci_lower, ci_upper = inf$ci_upper,
        n_bootstrap_failed = B - inf$n_bootstrap_used,
        ## Unadjusted contrasts, for reading the adjusted estimate against.
        mean_annotated = if (t$binary) mean(yy[x == 1]) else NA_real_,
        mean_other = if (t$binary) mean(yy[x == 0]) else NA_real_,
        frac_gain_with_age_annotated = if (t$binary) mean(g[x == 1] > 0) else NA_real_,
        frac_gain_with_age_other = if (t$binary) mean(g[x == 0] > 0) else NA_real_,
        covariates = paste(t$ax$covariates, collapse = ","))
}))
if (nrow(ann_res)) {
    ann_res[, q := stats::p.adjust(p, "BH"), by = .(outcome, adjustment)]
    ann_res[, estimate_meaning := fifelse(outcome == "debiased_sq_beta",
        "adjusted difference in debiased squared age slope, in units of the region mean (per SD if continuous)",
        "adjusted difference in signed age slope, methylation fraction per decade (per SD if continuous)")]
}
ann_res <- rbindlist(list(ann_res, rbindlist(ann_skipped)), fill = TRUE)
ann_res[, `:=`(cohort = cohort, region = region, run_id = opts$run_id, spec = "primary",
               role = "descriptive", n_bootstrap = B,
               inference = "donor_bootstrap_var_plus_chromosome_jackknife_var")]

write_atomic(res, file.path(run_dir, "results", "axis-tests.tsv"))
write_atomic(skipped, file.path(run_dir, "results", "axis-arms-skipped.tsv"))
write_atomic(ann_res, file.path(run_dir, "results", "annotation-age-associations.tsv"))
saveRDS(boot_ann, file.path(run_dir, "checkpoint", "annotation-bootstrap.rds"))
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
    axis_arms_not_fitted = paste(skipped$arm, collapse = ","),
    n_annotation_tests = as.character(length(ann_tests)),
    n_bootstrap_completed = as.character(B),
    max_bootstrap_failed = as.character(max(res$n_bootstrap_failed))
))
print(res[, .(spec, arm, outcome, n_vmrs, estimate = signif(estimate, 3),
              se = signif(se, 3), se_boot = signif(se_bootstrap, 3),
              se_jk = signif(se_jackknife, 3), p = signif(p, 3))])
