#!/usr/bin/env Rscript
#### 02b_greml_simulation_benchmark -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --arm ar1 [--allow-unlocked]
##   Rscript _h/00_new_run.R --arm observed --cohort AA --region dlpfc [--allow-unlocked]
##
## Arm 1 (ar1) is one run over every simulated sample size:
##   greml-sim-ar1-{YYYYMMDD}
## Arm 2 (observed) is one run per region:
##   greml-AA-{region}-{YYYYMMDD}
##
## Writes the task manifests every array stage reads, and prints the run ID as
## its last line for the submit driver.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02b_greml_simulation_benchmark", "_h")),
                 "greml_functions.R"))

MODULE <- "02b_greml_simulation_benchmark"

args <- commandArgs(trailingOnly = TRUE)
opts <- parse_v2_args(args, require = "arm")
allow_unlocked <- isTRUE(opts$allow_unlocked)
if (!opts$arm %in% c("ar1", "observed")) stop("--arm must be ar1 or observed")

cfg <- load_config("greml_benchmark")
assert_locked(list(greml_benchmark = cfg), allow_unlocked = allow_unlocked)
assert_greml_interpretation(cfg)
gcta_v <- assert_gcta_pin(cfg)
thr <- load_config("thresholds")

if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked); ",
         "production run IDs are derived, not chosen.")
}
module_root <- file.path(repo_root(), MODULE)
modes <- rbindlist(lapply(cfg$reml_modes, function(m)
    data.table(reml_mode = m$id, role = m$role, constrained = isTRUE(m$constrained),
               args = paste(unlist(m$args), collapse = " "))))
if (sum(modes$role == "primary") != 1L) stop("Exactly one reml_mode must be primary")

common_extra <- list(
    smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
    config_greml_benchmark_sha256 = file_sha256(file.path(repo_root(), "config",
                                                          "greml_benchmark.yml")),
    gcta_binary = cfg$gcta$binary,
    gcta_version = gcta_v,
    reml_primary_mode = modes[role == "primary", reml_mode],
    seed_namespace = cfg$seeds$namespace,
    simulated_phenotypes_only = "TRUE",
    absolute_pve_interpretation_allowed_for_observed_loci = "FALSE",
    reads_module_02 = "FALSE"
)

## ===================================================================== arm 1
if (opts$arm == "ar1") {
    a <- cfg$arm1
    sizes <- if (allow_unlocked) unlist(a$smoke$sample_sizes) else unlist(a$sample_sizes)
    n_pheno <- if (allow_unlocked) as.integer(a$smoke$n_phenotypes) else
        as.integer(a$n_phenotypes)
    root_in <- file.path(repo_root(), a$inputs_root)
    for (N in sizes) {
        d <- file.path(root_in, sprintf("sim_%d_indiv", N))
        need <- c(file.path(d, "plink_sim", paste0("simulated.", c("bed", "bim", "fam"))),
                  file.path(d, "simulated.phen"), file.path(d, "snp_phenotype_mapping.tsv"))
        miss <- need[!file.exists(need)]
        if (length(miss)) stop("Missing arm 1 input(s): ", paste(miss, collapse = ", "))
    }
    run <- new_run(module = cfg$module_tag, cohort = a$run_cohort,
                   region = a$run_region, module_root = module_root,
                   run_id = opts$run_id,
                   extra = c(common_extra, list(
                       arm = a$id,
                       sample_sizes = paste(sizes, collapse = ","),
                       n_phenotypes = as.character(n_pheno),
                       ld_score_region_kb = as.character(a$ld_score_region_kb),
                       n_ld_strata = as.character(a$n_ld_strata),
                       inputs_root = a$inputs_root)))
    for (s in c("config", "inputs", "work", "results/ldscore", "results/stratify",
                "results/reml", "status", "summary", "figures")) {
        dir.create(file.path(run$dir, s), recursive = TRUE, showWarnings = FALSE)
    }
    ld <- CJ(N = as.integer(sizes), chr = seq_len(as.integer(a$chromosomes)))
    ld[, task := .I]
    strat <- data.table(task = seq_along(sizes), N = as.integer(sizes))
    chunks <- rbindlist(lapply(sizes, function(N) {
        from <- seq(1L, n_pheno, by = as.integer(a$phenotypes_per_task))
        data.table(N = as.integer(N), pheno_from = from,
                   pheno_to = pmin(from + as.integer(a$phenotypes_per_task) - 1L, n_pheno))
    }))
    chunks[, task := .I]
    write_atomic(ld, file.path(run$dir, "config", "ldscore-tasks.tsv"))
    write_atomic(strat, file.path(run$dir, "config", "stratify-tasks.tsv"))
    write_atomic(chunks, file.path(run$dir, "config", "reml-tasks.tsv"))
    write_atomic(modes, file.path(run$dir, "config", "reml-modes.tsv"))
    append_manifest(run, list(n_ldscore_tasks = nrow(ld),
                              n_stratify_tasks = nrow(strat),
                              n_reml_tasks = nrow(chunks),
                              n_expected_reml_units = length(sizes) * n_pheno * nrow(modes)))
    cat(run$run_id, "\n", sep = "")
    quit(save = "no")
}

## ===================================================================== arm 2
a <- cfg$arm2
if (is.null(opts$cohort) || is.null(opts$region)) {
    stop("--arm observed needs --cohort and --region")
}
if (!identical(opts$cohort, a$cohort)) stop("arm 2 is prespecified for cohort ", a$cohort)
if (!opts$region %in% unlist(a$regions)) stop("region not in config arm2.regions")

up <- require_accepted_upstream(a$upstream, opts$cohort, opts$region,
                                allow_unaccepted = allow_unlocked)
if (is.na(up$run_id %||% NA_character_)) {
    stop("Module 01 has no accepted run for ", opts$cohort, " x ", opts$region,
         ". Arm 2 reads the run directory by ID, so it has no fallback.")
}
vmr_run_dir <- file.path(repo_root(), a$upstream, "_m", "runs", up$run_id)
design_n <- load_config("cohorts")$donor_counts[[opts$cohort]][[opts$region]]$design_n
if (is.null(design_n)) stop("No locked design_n for ", opts$cohort, " x ", opts$region,
                            " in config/cohorts.yml donor_counts")

catalog <- fread(file.path(vmr_run_dir, "vmr", "vmr_catalog.tsv"))
if (!identical(unique(catalog$vmr_set_id), up$vmr_set_id)) {
    stop("vmr_set_id in the catalog does not match the accepted row")
}
catalog <- catalog[!toupper(sub("^chr", "", chr)) %in% c("X", "Y")]
window_bp <- as.integer(thr$cis$window_bp)
min_cis <- as.integer(thr$cis$min_cis_variants)

## Pre-QC cis-window SNP count from each locus .bim. Cheap, deterministic, and
## needs no genotype read; the post-QC count is applied at task time.
bim_count <- function(chr, start, end) {
    lab <- sub("^chr", "", chr)
    bim <- file.path(vmr_run_dir, "plink_format", paste0("chr_", lab),
                     sprintf("TOPMed_LIBD-%s.%d_%d.bim", opts$cohort, start, end))
    if (!file.exists(bim)) return(0L)
    b <- fread(bim, header = FALSE, select = c(1L, 4L))
    as.integer(sum(sub("^chr", "", b[[1]]) == lab &
                   b[[2]] >= max(1L, start - window_bp) & b[[2]] <= end + window_bp))
}
catalog[, bim_snps := mapply(bim_count, chr, start, end)]
eligible <- catalog[bim_snps >= min_cis]
eligible[, stratum := cut(bim_snps,
                          breaks = unique(stats::quantile(bim_snps,
                                                          seq(0, 1, length.out = a$locus_strata + 1))),
                          include.lowest = TRUE, labels = FALSE)]
per <- if (allow_unlocked) as.integer(a$smoke$loci_per_stratum) else
    as.integer(a$loci_per_stratum)
set.seed(seed_for(cfg$seeds$namespace, opts$region, "locus_sample"))
loci <- eligible[, .SD[sample(.N, min(.N, per))], by = stratum]
setorder(loci, stratum, chr, start)
loci[, locus_index := .I]
loci[, task := ceiling(locus_index / as.integer(a$loci_per_task))]

reps <- if (allow_unlocked) as.integer(a$smoke$replicates_per_cell) else
    as.integer(a$replicates_per_cell)
grid <- CJ(h2 = as.numeric(unlist(a$h2_values)),
           architecture = as.character(unlist(a$architectures)),
           replicate = seq_len(reps))

run <- new_run(module = cfg$module_tag, cohort = opts$cohort, region = opts$region,
               module_root = module_root, run_id = opts$run_id,
               vmr_set_id = up$vmr_set_id,
               upstream = list(vmr_catalog_run_id = up$run_id),
               extra = c(common_extra, list(
                   arm = a$id,
                   design_n = as.character(design_n),
                   cis_window_bp = as.character(window_bp),
                   min_cis_variants = as.character(min_cis),
                   n_catalog_autosomal = as.character(nrow(catalog)),
                   n_eligible_pre_qc = as.character(nrow(eligible)),
                   n_loci_drawn = as.character(nrow(loci)),
                   replicates_per_cell = as.character(reps),
                   n_scenarios_per_locus = as.character(nrow(grid)))))
for (s in c("config", "work", "results/reml", "status", "summary", "figures")) {
    dir.create(file.path(run$dir, s), recursive = TRUE, showWarnings = FALSE)
}
write_atomic(loci[, .(task, locus_index, vmr_id, chrom = chr, start, end,
                      bim_snps, stratum)],
             file.path(run$dir, "config", "locus-tasks.tsv"))
write_atomic(grid, file.path(run$dir, "config", "scenario-grid.tsv"))
write_atomic(modes, file.path(run$dir, "config", "reml-modes.tsv"))
append_manifest(run, list(n_reml_tasks = max(loci$task),
                          n_expected_reml_units = nrow(loci) * nrow(grid) * nrow(modes)))
cat(run$run_id, "\n", sep = "")
