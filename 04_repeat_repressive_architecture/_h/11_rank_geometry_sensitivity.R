#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- rank geometry sensitivity ####
##
## Usage:
##   Rscript _h/11_rank_geometry_sensitivity.R --cohort AA
##   Rscript _h/11_rank_geometry_sensitivity.R --cohort AA --allow-unlocked --run-id <smoke id>
##
## A NON-GATING sensitivity over the three accepted cells, sealed as its own
## cross-region run (`rra-{cohort}-crossregion-{date}`). It asks whether the
## chromatin and repeat associations of the primary model could come from the
## locus geometry that the frozen joint model of Module 02 enters into the
## rank, rather than from the estimator outputs.
##
## Why it exists. The joint model adds log num_snps, log p_eff and ld_metric
## LINEARLY to its estimate (00_shared/rank_geometry.R). The primary model here
## adjusts for tested_snp_count but not for p_eff or LD. If effective rank or LD
## differed between chromatin states, the rank could inherit a chromatin
## association from those terms alone. A 2026-10-12 scratch refit found no such
## effect, but it was unsealed and not citable; this stage makes it a run.
##
## Two refits of the primary model, each beside the reproduced primary:
##   geometry_free_rank  the predictor is the rank of the joint estimate with
##                       its geometry part removed (computed per cell by
##                       00_shared/rank_geometry.R from the accepted Module 02
##                       run that the cell itself consumed);
##   adjust_geometry     the primary predictor, with the joint model's three
##                       within-cell geometry terms -- log(num_snps),
##                       log(p_eff), ld_metric -- added to the adjustment set.
##
## The model is the primary model of _h/02_test_association.R, rebuilt from the
## same config as in _h/09_shared_unique_split.R and checked, before any refit,
## to reproduce every sealed primary estimate to 1e-6. If it does not, the stage
## stops.
##
## What it is NOT. It changes no gate, claim, q-value or decision token. Every
## row carries gating = FALSE. Caudate LINE/L1 stays set aside as in the
## primary. The per-cell ATAC breakdowns are not refitted; the ATAC union is.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(repo_root(), "00_shared", "rank_geometry.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(sandwich)
    library(lmtest)
})

MODULE <- "04_repeat_repressive_architecture"
MODULE_TAG <- "rra"
LGV_MODULE <- "02_local_genetic_variance"
PRED <- "local_snp_contribution_score_z"
FREE <- "geometry_free_rank_z"
GEOMETRY_ADJUST <- c("log(num_snps)", "log(p_eff)", "ld_metric")
REPRO_TOL <- 1e-6
JOINT_SETTINGS <- file.path(repo_root(), LGV_MODULE, "config", "joint-pve-20260820.tsv")

opts <- parse_v2_args(require = "cohort")
allow_unlocked <- isTRUE(opts$allow_unlocked)
cohort <- opts$cohort
if (!identical(parse_cell(cohort)$cell_kind, "arm")) {
    stop("The sensitivity runs on a discovery arm's cells, not the cell '", cohort, "'")
}
if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked)")
}

annot <- load_config("repeat_annotations")
lgc_cfg <- load_config("local_genetic_control")
assert_locked(list(repeat_annotations = annot, local_genetic_control = lgc_cfg),
              allow_unlocked = allow_unlocked)
settings <- read_joint_settings_table(JOINT_SETTINGS)
regions <- c("caudate", "dlpfc", "hippocampus")

## ------------------------------------------------------------ the three cells
up <- lapply(regions, function(re)
    require_accepted_upstream(MODULE, cohort, re, allow_unaccepted = FALSE))
names(up) <- regions
cell_dir <- function(re) file.path(repo_root(), MODULE, "_m", "runs", up[[re]]$run_id)
cell_manifest <- function(re, field) {
    m <- fread(file.path(cell_dir(re), "manifest.tsv"), colClasses = "character")
    v <- m$value[m$field == field]
    if (length(v) != 1L) stop(up[[re]]$run_id, " manifest lacks a unique ", field)
    v
}

## The config the cells were fitted under must say what the config this stage
## reads says. A byte checksum is too strict for that: commit 912de5d45 edited
## only comments in config/repeat_annotations.yml after the cells were sealed.
## So: recover the version of the file the cell checksummed, and require its
## PARSED content to equal the current file's. A comment-only edit passes; any
## change to a setting stops the stage. The reproduction guard below then
## checks the fitted model itself.
##
## The version is found by checksum in the file's history, not at the cell's
## recorded git_commit: the history was rewritten after the 2026-09-25 cells
## were sealed, so their recorded commit (d61f83b4c) no longer exists here.
cfg_sha <- attr(annot, "config_sha256")
cfg_rel <- "config/repeat_annotations.yml"
config_version_by_sha <- function(sha) {
    commits <- system2("git", c("-C", repo_root(), "log", "--format=%H", "--", cfg_rel),
                       stdout = TRUE)
    for (commit in commits) {
        tmp <- tempfile(fileext = ".yml")
        status <- system2("git", c("-C", repo_root(), "show", paste0(commit, ":", cfg_rel)),
                          stdout = tmp, stderr = FALSE)
        if (identical(as.integer(status), 0L) && identical(file_sha256(tmp), sha)) {
            return(list(commit = commit, path = tmp))
        }
    }
    NULL
}
config_provenance <- list()
for (re in regions) {
    rec <- cell_manifest(re, "config_repeat_annotations_sha256")
    if (identical(rec, cfg_sha)) {
        config_provenance[[re]] <- "identical"
        next
    }
    hit <- config_version_by_sha(rec)
    if (is.null(hit)) {
        stop(up[[re]]$run_id, " was fitted under config_repeat_annotations_sha256 ",
             rec, ", which is not any version of ", cfg_rel, " in this history.")
    }
    if (!identical(yaml::read_yaml(hit$path), yaml::read_yaml(attr(annot, "config_file")))) {
        stop(up[[re]]$run_id, " was fitted under a config whose settings differ from ",
             "the current ", cfg_rel, ". Run this stage from commit ", hit$commit, ".")
    }
    config_provenance[[re]] <- paste0("settings_identical_to_", substr(hit$commit, 1, 9))
    message("[04] ", re, ": config differs from the cell's version (",
            substr(hit$commit, 1, 9), ") only outside its settings")
}

## The Module 02 run each cell consumed must still be the accepted one, so the
## geometry-free rank is built from the same estimates the cell ranked on.
lgv_run <- vapply(regions, function(re) {
    rid <- cell_manifest(re, "upstream_local_genetic_variance_run_id")
    acc <- require_accepted_upstream(LGV_MODULE, cohort, re, allow_unaccepted = FALSE)$run_id
    if (!identical(rid, acc)) {
        stop(up[[re]]$run_id, " consumed ", rid, " but the accepted Module 02 run is ", acc)
    }
    rid
}, character(1))

feats <- lapply(regions, function(re) {
    f <- fread(file.path(cell_dir(re), "results", "vmr-features.tsv"))
    lgc <- load_local_genetic_control(lgv_run[[re]], re, cohort, eligible_only = TRUE)
    shas <- unique(as.character(lgc$joint_model_sha256))
    if (!identical(shas, lgc_cfg$validation$frozen_model_sha256)) {
        stop(lgv_run[[re]], " was not scored by the frozen joint model")
    }
    dec <- decompose_rank_geometry(lgc, settings, cell = cohort, region = re)
    geo <- merge(dec$per_vmr[, .(vmr_id, z_lgv = local_snp_contribution_score_z,
                                 geometry_free_rank_z)],
                 lgc[, .(vmr_id, num_snps, p_eff, ld_metric)], by = "vmr_id")
    out <- merge(f, geo, by = "vmr_id")
    if (nrow(out) != nrow(f)) {
        stop(re, ": ", nrow(f) - nrow(out), " feature rows have no Module 02 geometry")
    }
    z_err <- max(abs(out[[PRED]] - out$z_lgv))
    if (!is.finite(z_err) || z_err > 1e-10) {
        stop(re, ": the cell's rank differs from ", lgv_run[[re]], " (", z_err, ")")
    }
    out[, z_lgv := NULL][]
})
names(feats) <- regions
sealed <- lapply(regions, function(re)
    fread(file.path(cell_dir(re), "results", "association-results.tsv"))[
        analysis_set == "primary" & predictor == PRED])
names(sealed) <- regions

## ---------------------------------------------------------- the primary model
## Mirrors _h/02_test_association.R and _h/09_shared_unique_split.R.
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
role_of <- function(o) {
    if (o %in% BH_FAMILY) return("bh_family")
    annot$multiple_testing$outside_family[[o]]$role %||% "control"
}
## The family and every control except the per-cell-type ATAC breakdowns, which
## are a secondary split of atac_union_frac and are refitted through it.
OUTCOMES <- Filter(function(o) role_of(o) != "celltype_breakdown_secondary",
                   c(BH_FAMILY, CONTROLS))
set_aside <- annot$interpretation$technically_confounded_regions

fit_hc3 <- function(d, outcome, predictor, terms) {
    form <- stats::as.formula(paste(outcome, "~", predictor, "+",
                                    paste(terms, collapse = " + ")))
    fam <- stats::quasibinomial(link = annot$primary_model$link %||% "logit")
    fit <- tryCatch(suppressWarnings(stats::glm(form, data = d, family = fam)),
                    error = function(e) NULL)
    if (is.null(fit) || is.na(stats::coef(fit)[predictor])) {
        return(list(estimate = NA_real_, se = NA_real_, p = NA_real_, n_fitted = 0L))
    }
    ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(
        fit, type = annot$primary_model$vcov %||% "HC3"))
    list(estimate = ct[predictor, 1], se = ct[predictor, 2], p = ct[predictor, 4],
         n_fitted = as.integer(stats::nobs(fit)))
}

## ------------------------------------------------------- reproduction guard
repro <- rbindlist(lapply(regions, function(re) {
    terms <- covariates_for(re)
    rbindlist(lapply(OUTCOMES, function(o) {
        r <- fit_hc3(feats[[re]], o, PRED, terms)
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
    stop(nrow(bad), " primary fit(s) do not reproduce the sealed cells.")
}
message("[04] reproduction guard: ", nrow(repro), " primary fits reproduce to <= ",
        REPRO_TOL)

## ------------------------------------------------------------- the refits
MODELS <- list(
    primary = function(terms) list(predictor = PRED, terms = terms),
    geometry_free_rank = function(terms) list(predictor = FREE, terms = terms),
    adjust_geometry = function(terms) list(predictor = PRED,
                                           terms = c(terms, GEOMETRY_ADJUST)))
assoc <- rbindlist(lapply(regions, function(re) {
    base_terms <- covariates_for(re)
    rbindlist(lapply(OUTCOMES, function(o) {
        rbindlist(lapply(names(MODELS), function(mname) {
            spec <- MODELS[[mname]](base_terms)
            r <- fit_hc3(feats[[re]], o, spec$predictor, spec$terms)
            data.table(region = re, outcome = o, model = mname,
                       predictor = spec$predictor,
                       added_terms = if (mname == "adjust_geometry")
                           paste(GEOMETRY_ADJUST, collapse = " + ") else NA_character_,
                       n_fitted = r$n_fitted, estimate = r$estimate,
                       se_hc3 = r$se, p_hc3 = r$p)
        }))
    }))
}))
prim <- assoc[model == "primary", .(region, outcome, estimate_primary = estimate)]
assoc <- merge(assoc, prim, by = c("region", "outcome"))
assoc[, `:=`(ratio_to_primary = estimate / estimate_primary,
             same_sign_as_primary = sign(estimate) == sign(estimate_primary))]
assoc[, outcome_role := vapply(outcome, role_of, character(1))]
assoc[, set_aside_from_claim := mapply(function(o, re) re %in% unlist(set_aside[[o]]),
                                       outcome, region)]
assoc[, `:=`(gating = FALSE, sensitivity = "rank_geometry",
             population = cohort,
             source_run_id = vapply(region, function(re) up[[re]]$run_id, character(1)),
             score_run_id = unname(lgv_run[region]))]
setorder(assoc, region, outcome, model)

## ------------------------------------------------------------------ the run
run <- new_run(
    module = MODULE_TAG, cohort = cohort, region = "crossregion",
    module_root = file.path(repo_root(), MODULE), run_id = opts$run_id,
    upstream = c(stats::setNames(lapply(regions, function(re) up[[re]]$run_id),
                                 paste0("repeat_architecture_", regions)),
                 stats::setNames(as.list(lgv_run),
                                 paste0("local_genetic_variance_", regions))),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_repeat_annotations_sha256 = cfg_sha,
        config_repeat_annotations_vs_cells = paste(
            paste0(regions, "=", unlist(config_provenance[regions])), collapse = ";"),
        config_local_genetic_control_sha256 = attr(lgc_cfg, "config_sha256"),
        joint_settings_sha256 = file_sha256(JOINT_SETTINGS),
        sensitivity = "rank_geometry",
        gating = "FALSE",
        geometry_adjustment_terms = paste(GEOMETRY_ADJUST, collapse = " + "),
        reproduction_tolerance = REPRO_TOL))
res_dir <- file.path(run$dir, "results")
dir.create(res_dir, showWarnings = FALSE)
code_dir <- file.path(run$dir, "code")
dir.create(code_dir, showWarnings = FALSE)
invisible(file.copy(c(file.path(repo_root(), MODULE, "_h", "11_rank_geometry_sensitivity.R"),
                      file.path(repo_root(), "00_shared", "rank_geometry.R")), code_dir))

write_atomic(repro, file.path(res_dir, "geometry-reproduction-check.tsv"))
write_atomic(assoc, file.path(res_dir, "geometry-sensitivity-association.tsv"))
writeLines(c(
    "Interpretation constraints carried by this run:",
    "  - NON-GATING. No gate, claim, q-value or decision token of the accepted",
    "    cells changes.",
    "  - geometry_free_rank is a sensitivity predictor built from the frozen",
    "    joint model's estimator part; it is not a replacement endpoint.",
    "  - Agreement with the primary shows the associations are not produced by",
    "    the joint model's geometry terms. It does not exclude other dependence",
    "    of the estimators themselves on local genotype structure.",
    "  - Caudate LINE/L1 remains set aside from the claim (set_aside_from_claim)."
), file.path(res_dir, "interpretation-constraints.txt"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "session-info.txt"))

if (!allow_unlocked) close_run(run)
message("[04] geometry run ", run$run_id, if (allow_unlocked) " (smoke, not sealed)" else " sealed")
print(dcast(assoc, region + outcome ~ model, value.var = "estimate"), digits = 3)
cat(run$run_id, "\n", sep = "")
