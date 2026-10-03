#!/usr/bin/env Rscript
#### 02b arm 2 -- cis-window REML on real genotypes, simulated phenotypes ####
##
## For each locus in the task: load the real cis-window genotype and the locked
## covariates through 00_shared/locus_io.R::load_observed_locus() (same window,
## SNP QC, donor alignment and expected n as every v2 consumer), write its GRM,
## then for every scenario simulate a phenotype at a known h2 with
## simulate_phenotype_on_observed_genotype() and fit GCTA REML in each REML
## mode. The observed methylation phenotype is never fitted.
##
## Usage: Rscript _h/05_arm2_greml.R --run-id <id> [--task k]

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02b_greml_simulation_benchmark", "_h")),
                 "greml_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
cfg <- load_greml_config(run_dir)
man <- read_run_manifest(run_dir)
task_id <- as.integer(opts$task %||% Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(task_id)) stop("No --task and no SLURM_ARRAY_TASK_ID")
threads <- as.integer(Sys.getenv("V2_CPUS", "1"))

cohort <- man[["cohort"]]; region <- man[["region"]]
vmr_run_dir <- file.path(repo_root(), cfg$arm2$upstream, "_m", "runs",
                         man[["upstream_vmr_catalog_run_id"]])
loci <- fread(file.path(run_dir, "config", "locus-tasks.tsv"))[task == task_id]
if (!nrow(loci)) stop("No loci for task ", task_id)
grid <- fread(file.path(run_dir, "config", "scenario-grid.tsv"))
modes <- fread(file.path(run_dir, "config", "reml-modes.tsv"))
patterns <- unlist(cfg$gcta$estimation_failure_patterns)
tmp <- file.path(run_dir, "work", sprintf("task_%d", task_id))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

unit_id <- function(vmr, g, m) arm2_unit_id(vmr, g$h2, g$architecture,
                                           g$replicate, m)
rows <- list(); st <- list(); geo <- list()
for (i in seq_len(nrow(loci))) {
    L <- loci[i]
    task <- list(chrom = L$chrom, start = L$start, end = L$end, vmr_id = L$vmr_id)
    loc <- load_observed_locus(task, cohort = cohort, vmr_run_dir = vmr_run_dir,
                               min_cis_variants = as.integer(man[["min_cis_variants"]]),
                               expected_n = as.integer(man[["design_n"]]),
                               backing_tag = "greml")
    if (!identical(loc$status, "ok")) {
        for (j in seq_len(nrow(grid))) for (m in modes$reml_mode) {
            st[[length(st) + 1L]] <- data.table(
                unit = unit_id(L$vmr_id, grid[j], m), status = "qc_failed",
                reason = paste0("locus_", loc$status, ": ", loc$reason))
        }
        geo[[length(geo) + 1L]] <- data.table(vmr_id = L$vmr_id, locus_status = loc$status,
                                              locus_reason = loc$reason, n_donors = NA_integer_,
                                              n_snps_qc = NA_integer_, bim_snps = L$bim_snps)
        next
    }
    G <- loc$genotype
    grm <- grm_from_dosage(G)
    ids <- loc$metadata[, c("FID", "IID")]
    pre <- file.path(tmp, "locus")
    write_grm(grm, ids, pre)
    qcov <- file.path(tmp, "qcovar.txt")
    fwrite(cbind(as.data.table(ids), as.data.table(loc$covariates)), qcov,
           sep = "\t", col.names = FALSE)
    geo[[length(geo) + 1L]] <- data.table(vmr_id = L$vmr_id, locus_status = "ok",
                                          locus_reason = NA_character_,
                                          n_donors = nrow(G), n_snps_qc = grm$n_snps,
                                          bim_snps = L$bim_snps)
    for (j in seq_len(nrow(grid))) {
        g <- grid[j]
        set.seed(seed_for(cfg$seeds$namespace, region, L$vmr_id,
                          paste(as.character(g$h2), g$architecture, g$replicate)))
        sim <- simulate_phenotype_on_observed_genotype(G, loc$covariates, g$h2,
                                                       g$architecture)
        phen <- file.path(tmp, "pheno.txt")
        fwrite(cbind(as.data.table(ids), y = sim$phenotype), phen,
               sep = "\t", col.names = FALSE)
        for (k in seq_len(nrow(modes))) {
            m <- modes[k]
            out <- file.path(tmp, paste0("fit_", m$reml_mode))
            unlink(paste0(out, c(".hsq", ".gcta.log", ".log")))
            margs <- if (nzchar(m$args)) strsplit(m$args, " ", fixed = TRUE)[[1]] else character(0)
            r <- run_gcta(cfg$gcta$binary,
                          c("--reml", "--grm", pre, "--pheno", phen, "--qcovar", qcov,
                            "--reml-maxit", cfg$gcta$reml_maxit, margs),
                          out, threads = threads)
            cl <- classify_reml(r$exit, paste0(out, ".hsq"), r$log, patterns,
                                reml_divergence_rule(cfg))
            u <- unit_id(L$vmr_id, g, m$reml_mode)
            st[[length(st) + 1L]] <- data.table(unit = u, status = cl$status, reason = cl$reason)
            if (cl$status == "completed") {
                rows[[length(rows) + 1L]] <- cbind(
                    data.table(arm = man[["arm"]], cohort = cohort, region = region,
                               vmr_id = L$vmr_id, stratum = L$stratum,
                               n_donors = nrow(G), n_snps_qc = grm$n_snps,
                               h2_nominal = g$h2, architecture = g$architecture,
                               replicate = g$replicate, realized_h2 = sim$realized_h2,
                               reml_mode = m$reml_mode, constrained = m$constrained),
                    parse_hsq(paste0(out, ".hsq")))
            } else if (cl$status == "failed") {
                file.copy(r$log, file.path(run_dir, "logs",
                                           sprintf("failed_task%d_%d.gcta.log", task_id,
                                                   length(st))))
            }
        }
    }
    unlink(list.files(tmp, full.names = TRUE))
}
write_atomic(rbindlist(rows, fill = TRUE),
             file.path(run_dir, "results", "reml", sprintf("task_%d.tsv", task_id)))
write_atomic(rbindlist(geo, fill = TRUE),
             file.path(run_dir, "results", "reml", sprintf("loci_task_%d.tsv", task_id)))
st <- rbindlist(st)
write_status(st, file.path(run_dir, "status", sprintf("reml_task_%d.tsv", task_id)))
unlink(tmp, recursive = TRUE)
if (nrow(st[status == "failed"])) {
    stop(nrow(st[status == "failed"]), " REML call(s) failed computationally; first: ",
         st[status == "failed"]$reason[1])
}
