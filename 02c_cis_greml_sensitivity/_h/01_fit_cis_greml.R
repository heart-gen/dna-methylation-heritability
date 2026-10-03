#!/usr/bin/env Rscript
#### 02c -- cis-window GREML on the observed VMR phenotypes, one task chunk ####
##
## For each VMR: the observed phenotype, genotype and covariates come from
## 00_shared/locus_io.R::load_observed_locus() -- the reader Module 02 Stage 01
## uses -- so the cis window, SNP QC, donor alignment, design n and covariate
## matrix are Module 02's by construction. The QC'd SNP count is then CHECKED
## against Module 02's own num_snps for the VMR, so "the same SNPs" is a
## measured property of the run, not an assumption. One cis-window GRM is
## written from exactly those SNPs (00_shared/gcta.R::grm_from_dosage, which
## matches gcta --make-grm) and REML is fitted in every configured mode with the
## covariates as --qcovar. Estimate, SE, logL, LRT and convergence are kept;
## nothing is classified on them.
##
## Usage: Rscript _h/01_fit_cis_greml.R --run-id <id> [--task k]

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02c_cis_greml_sensitivity", "_h")),
                 "cgs_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
cfg <- load_cgs_config(run_dir)
man <- read_run_manifest(run_dir)
task_id <- as.integer(opts$task %||% Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(task_id)) stop("No --task and no SLURM_ARRAY_TASK_ID")
threads <- as.integer(Sys.getenv("V2_CPUS", "1"))

cohort <- man[["cohort"]]
vmr_run_dir <- file.path(repo_root(), cfg$upstreams$vmr_catalog, "_m", "runs",
                         man[["upstream_vmr_catalog_run_id"]])
tasks <- fread(file.path(run_dir, "config", "vmr-tasks.tsv"))[task == task_id]
if (!nrow(tasks)) stop("No VMRs for task ", task_id)
modes <- fread(file.path(run_dir, "config", "reml-modes.tsv"))
patterns <- unlist(cfg$gcta$estimation_failure_patterns)
divergence <- reml_divergence_rule(cfg)
tmp <- file.path(run_dir, "work", sprintf("task_%d", task_id))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

rows <- list(); st <- list(); loci <- list()
for (i in seq_len(nrow(tasks))) {
    V <- tasks[i]
    loc <- load_observed_locus(list(chrom = V$chrom, start = V$start, end = V$end,
                                    vmr_id = V$vmr_id),
                               cohort = cohort, vmr_run_dir = vmr_run_dir,
                               min_cis_variants = as.integer(man[["min_cis_variants"]]),
                               expected_n = as.integer(man[["design_n"]]),
                               backing_tag = "cgs")
    if (!identical(loc$status, "ok")) {
        for (m in modes$reml_mode) st[[length(st) + 1L]] <- data.table(
            unit = unit_id(V$vmr_id, m), status = "qc_failed",
            reason = paste0("locus_", loc$status, ": ", loc$reason))
        loci[[length(loci) + 1L]] <- data.table(vmr_id = V$vmr_id, locus_status = loc$status,
                                                locus_reason = loc$reason, n_donors = NA_integer_,
                                                n_snps_qc = NA_integer_,
                                                module02_num_snps = V$module02_num_snps,
                                                snp_count_matches_module02 = NA)
        next
    }
    n_snps <- ncol(loc$genotype)
    match02 <- if (is.finite(V$module02_num_snps)) n_snps == V$module02_num_snps else NA
    loci[[length(loci) + 1L]] <- data.table(vmr_id = V$vmr_id, locus_status = "ok",
                                            locus_reason = NA_character_,
                                            n_donors = nrow(loc$genotype), n_snps_qc = n_snps,
                                            module02_num_snps = V$module02_num_snps,
                                            snp_count_matches_module02 = match02)
    grm <- grm_from_dosage(loc$genotype)
    ids <- loc$metadata[, c("FID", "IID")]
    pre <- file.path(tmp, "locus")
    write_grm(grm, ids, pre)
    qcov <- file.path(tmp, "qcovar.txt")
    fwrite(cbind(as.data.table(ids), as.data.table(loc$covariates)), qcov,
           sep = "\t", col.names = FALSE)
    phen <- file.path(tmp, "pheno.txt")
    fwrite(cbind(as.data.table(ids), y = loc$y), phen, sep = "\t", col.names = FALSE)
    for (k in seq_len(nrow(modes))) {
        m <- modes[k]
        out <- file.path(tmp, paste0("fit_", m$reml_mode))
        unlink(paste0(out, c(".hsq", ".gcta.log", ".log")))
        margs <- if (nzchar(m$args)) strsplit(m$args, " ", fixed = TRUE)[[1]] else character(0)
        r <- run_gcta(cfg$gcta$binary,
                      c("--reml", "--grm", pre, "--pheno", phen, "--qcovar", qcov,
                        "--reml-maxit", cfg$gcta$reml_maxit, margs),
                      out, threads = threads)
        cl <- classify_reml(r$exit, paste0(out, ".hsq"), r$log, patterns, divergence)
        st[[length(st) + 1L]] <- data.table(unit = unit_id(V$vmr_id, m$reml_mode),
                                            status = cl$status, reason = cl$reason)
        fit <- if (cl$status == "completed") parse_hsq(paste0(out, ".hsq")) else
            data.table(h2_hat = NA_real_, h2_se = NA_real_, vg_hat = NA_real_,
                       ve_hat = NA_real_, vp_hat = NA_real_, logL = NA_real_,
                       logL0 = NA_real_, lrt = NA_real_, pval = NA_real_, n_used = NA_real_)
        rows[[length(rows) + 1L]] <- cbind(
            data.table(vmr_id = V$vmr_id, chrom = V$chrom, start = V$start, end = V$end,
                       reml_mode = m$reml_mode, reml_role = m$role,
                       converged = cl$status == "completed", fit_status = cl$status,
                       fit_reason = cl$reason, n_donors = nrow(loc$genotype),
                       n_snps_qc = n_snps), fit)
        if (cl$status == "failed") file.copy(r$log, file.path(run_dir, "logs",
            sprintf("failed_task%d_%s_%s.gcta.log", task_id, gsub("[:]", "_", V$vmr_id),
                    m$reml_mode)))
    }
    unlink(list.files(tmp, full.names = TRUE))
}
write_atomic(rbindlist(rows, fill = TRUE),
             file.path(run_dir, "results", sprintf("fits_task_%d.tsv", task_id)))
write_atomic(rbindlist(loci, fill = TRUE),
             file.path(run_dir, "results", sprintf("loci_task_%d.tsv", task_id)))
st <- rbindlist(st)
write_status(st, file.path(run_dir, "status", sprintf("task_%d.tsv", task_id)))
unlink(tmp, recursive = TRUE)
if (nrow(st[status == "failed"])) {
    stop(nrow(st[status == "failed"]), " REML call(s) failed computationally; first: ",
         st[status == "failed"]$reason[1])
}
