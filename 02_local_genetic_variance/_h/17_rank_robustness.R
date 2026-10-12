#!/usr/bin/env Rscript
#### 02_local_genetic_variance -- Stage 17: rank robustness to geometry and genotype PCs ####
##
## Usage:
##   Rscript _h/17_rank_robustness.R
##   Rscript _h/17_rank_robustness.R --allow-unlocked --run-id <smoke id>
##
## A NON-GATING sensitivity over the accepted cells, sealed as one run
## (`lgv-arms-crossregion-{date}`). It reads sealed per-VMR tables only and
## re-estimates nothing. It asks two questions of the relative rank.
##
## 1. LOCUS GEOMETRY. The frozen joint model enters log num_snps, log p_eff and
##    ld_metric linearly (see 00_shared/rank_geometry.R). Do those terms, rather
##    than the estimator outputs, order the VMRs? For each of the six accepted
##    arm cells, the sealed estimate is split exactly into an estimator part and
##    a geometry part. The run reports the geometry part's share of variance and
##    the rank correlation between the sealed rank and a rank built from the
##    estimator part alone, and it writes that geometry-free rank per VMR.
##    04_repeat_repressive_architecture/_h/11_rank_geometry_sensitivity.R refits
##    the chromatin models on the same quantity, computed by the same function.
##
## 2. GENOTYPE PCs. The rank phenotype of the arm cells is adjusted for age, sex
##    and diagnosis only (00_shared/locus_io.R). It is NOT adjusted for snpPCs:
##    the snpPC residualization in 01_vmr_catalog feeds the CpG variability
##    cutoff, not the per-VMR phenotype, which 01_vmr_catalog/_h/02_summarize.R
##    writes from the smoothed methylation directly. The donor-group cell
##    all_individuals.AA scores the same Black American donors (the ordered
##    donor checksums are asserted equal) with three within-group snpPCs added.
##    On VMRs with identical coordinates in both catalogs, the rank correlation
##    between the two cells bounds how much the omitted PCs move the ordering.
##    It is an upper bound on the PC effect, because the donor-group cell also
##    differs in two other ways: its genotypes passed the pooled catalog's
##    quality control (fewer cis SNPs), and its CpG coverage filter was applied
##    across all donors of the pooled catalog. Each run also carries its own
##    Monte Carlo noise (elastic-net folds, BSLMM chains).
##
## What it is NOT. It changes no score, no eligibility, no gate and no claim.
## Every row carries gating = FALSE. The geometry-free rank is a sensitivity
## predictor, not a replacement endpoint.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(repo_root(), "00_shared", "rank_geometry.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "02_local_genetic_variance"
MODULE_TAG <- "lgv"
ARMS <- c("AA", "all_individuals")
DONOR_GROUP_CELL <- "all_individuals.AA"
REGIONS <- c("caudate", "dlpfc", "hippocampus")
JOINT_SETTINGS <- file.path(repo_root(), MODULE, "config", "joint-pve-20260820.tsv")

opts <- parse_v2_args(require = character())
allow_unlocked <- isTRUE(opts$allow_unlocked)
if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked)")
}

lgc_cfg <- load_config("local_genetic_control")
assert_locked(list(local_genetic_control = lgc_cfg), allow_unlocked = allow_unlocked)
frozen_sha <- lgc_cfg$validation$frozen_model_sha256
settings <- read_joint_settings_table(JOINT_SETTINGS)

run_dir_of <- function(run_id) file.path(repo_root(), MODULE, "_m", "runs", run_id)
manifest_value <- function(run_id, field) {
    m <- fread(file.path(run_dir_of(run_id), "manifest.tsv"), colClasses = "character")
    v <- m$value[m$field == field]
    if (length(v) != 1L) stop(run_id, " manifest lacks a unique ", field)
    v
}
check_frozen_model <- function(dt, run_id) {
    shas <- unique(as.character(dt$joint_model_sha256))
    if (!identical(shas, frozen_sha)) {
        stop(run_id, " was not scored by the frozen model ", frozen_sha, ": ",
             paste(shas, collapse = ", "))
    }
}

## ------------------------------------------------------ 1. locus geometry
up <- list()
decomp <- list()
for (arm in ARMS) for (re in REGIONS) {
    u <- require_accepted_upstream(MODULE, arm, re, allow_unaccepted = FALSE)
    up[[paste(arm, re, sep = "_")]] <- u$run_id
    dt <- load_local_genetic_control(u$run_id, re, arm, eligible_only = TRUE)
    check_frozen_model(dt, u$run_id)
    decomp[[paste(arm, re)]] <- decompose_rank_geometry(dt, settings, cell = arm, region = re)
    message("[17] ", arm, " ", re, ": ", nrow(dt), " eligible VMRs decomposed")
}
geometry_summary <- rbindlist(lapply(decomp, `[[`, "summary"))
geometry_coefficients <- rbindlist(lapply(decomp, `[[`, "coefficients"))
geometry_free <- rbindlist(lapply(decomp, `[[`, "per_vmr"))
assert_frozen_geometry(geometry_coefficients)

## ---------------------------------------------------- 2. genotype PCs
pc_rows <- list()
for (re in REGIONS) {
    primary_run <- up[[paste("AA", re, sep = "_")]]
    group_run <- require_accepted_upstream(MODULE, DONOR_GROUP_CELL, re,
                                           allow_unaccepted = FALSE)$run_id
    up[[paste(DONOR_GROUP_CELL, re, sep = "_")]] <- group_run
    if (!identical(manifest_value(primary_run, "ordered_donor_checksum"),
                   manifest_value(group_run, "ordered_donor_checksum"))) {
        stop(re, ": ", primary_run, " and ", group_run, " do not score the same ",
             "donors in the same order; the comparison would not isolate the PCs.")
    }
    ## The donor-group cell's per-VMR table is tracked as Stage 14's combined
    ## output; refuse it unless it was built from the accepted run.
    combined_f <- file.path(repo_root(), MODULE, "_m", "combined",
                            sprintf("local-genetic-control-donor-group-%s.tsv", re))
    grp <- fread(combined_f)
    if (!identical(unique(grp$upstream_lgv_run_id_1), group_run)) {
        stop(combined_f, " was not built from the accepted ", group_run)
    }
    check_frozen_model(load_local_genetic_control(primary_run, re, "AA",
                                                  eligible_only = TRUE), primary_run)
    pri <- load_local_genetic_control(primary_run, re, "AA", eligible_only = FALSE)
    pri <- pri[local_genetic_control_eligible %in% c(TRUE, "TRUE")]
    grp <- grp[local_genetic_control_eligible_AA %in% c(TRUE, "TRUE")]
    m <- merge(pri, grp, by = c("chrom", "start", "end"), suffixes = c("", ".grp"))
    sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))
    pc_rows[[re]] <- data.table(
        region = re,
        primary_run_id = primary_run, donor_group_run_id = group_run,
        n_donors = as.integer(manifest_value(primary_run, "n_donors")),
        ordered_donor_checksum_identical = TRUE,
        primary_covariates = "age + sex + diagnosis",
        donor_group_covariates = "age + sex + diagnosis + 3 within-group snpPCs",
        n_primary_eligible = nrow(pri),
        n_identical_vmrs_both_eligible = nrow(m),
        median_snp_count_ratio_group_to_primary = stats::median(m$num_snps_AA / m$num_snps),
        spearman_rank_z = sp(m$local_snp_contribution_score_z,
                             m$local_snp_contribution_score_z_AA),
        spearman_he_h2 = sp(m$he_h2, m$he_h2_AA),
        spearman_bslmm_pve = sp(m$bslmm_pve, m$bslmm_pve_AA),
        spearman_rho2_oof = sp(m$rho2_oof, m$rho2_oof_AA),
        spearman_r2_oof = sp(m$r2_oof, m$r2_oof_AA))
    message("[17] ", re, ": PC comparison on ", nrow(m), " identical VMRs")
}
pc_sensitivity <- rbindlist(pc_rows)

label <- function(dt) dt[, `:=`(gating = FALSE, sensitivity = "rank_geometry_and_genotype_pcs")][]
geometry_summary <- label(geometry_summary)
geometry_coefficients <- label(geometry_coefficients)
pc_sensitivity <- label(pc_sensitivity)

## ------------------------------------------------------------------ the run
run <- new_run(
    module = MODULE_TAG, cohort = "arms", region = "crossregion",
    module_root = file.path(repo_root(), MODULE), run_id = opts$run_id,
    upstream = up,
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        sensitivity = "rank_geometry_and_genotype_pcs",
        gating = "FALSE",
        config_local_genetic_control_sha256 = attr(lgc_cfg, "config_sha256"),
        joint_settings_sha256 = file_sha256(JOINT_SETTINGS),
        frozen_model_sha256 = frozen_sha))
res_dir <- file.path(run$dir, "results")
dir.create(res_dir, showWarnings = FALSE)
code_dir <- file.path(run$dir, "code")
dir.create(code_dir, showWarnings = FALSE)
invisible(file.copy(c(file.path(repo_root(), MODULE, "_h", "17_rank_robustness.R"),
                      file.path(repo_root(), "00_shared", "rank_geometry.R")), code_dir))

write_atomic(geometry_summary, file.path(res_dir, "geometry-decomposition-summary.tsv"))
write_atomic(geometry_coefficients, file.path(res_dir, "geometry-decomposition-coefficients.tsv"))
write_atomic(geometry_free, file.path(res_dir, "geometry-free-rank.tsv"))
write_atomic(pc_sensitivity, file.path(res_dir, "genotype-pc-sensitivity.tsv"))
writeLines(c(
    "Interpretation constraints carried by this run:",
    "  - NON-GATING. No score, eligibility, gate or claim changes.",
    "  - The geometry-free rank is a sensitivity predictor built from the frozen",
    "    model's estimator part. It is not a replacement endpoint and carries no",
    "    calibration of its own.",
    "  - The genotype-PC comparison is an upper bound on the effect of the omitted",
    "    snpPCs: the donor-group cell also used the pooled catalog's genotype QC",
    "    and CpG coverage filter, and each run carries its own Monte Carlo noise.",
    "  - Rank correlations compare orderings within a cell, never score levels."
), file.path(res_dir, "interpretation-constraints.txt"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "session-info.txt"))

if (!allow_unlocked) close_run(run)
message("[17] run ", run$run_id, if (allow_unlocked) " (smoke, not sealed)" else " sealed")
print(geometry_summary[, .(cell, region, n_eligible, var_share_geometry,
                           spearman_rank_vs_geometry_free_rank)], digits = 3)
print(pc_sensitivity[, .(region, n_identical_vmrs_both_eligible, spearman_rank_z,
                         median_snp_count_ratio_group_to_primary)], digits = 3)
cat(run$run_id, "\n", sep = "")
