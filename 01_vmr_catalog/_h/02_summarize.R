#### 01_vmr_catalog / 02_summarize: call VMRs and compute region methylation ####
##
## Merges the legacy 03.extract_vmr.R and 04.cal_vmr.R. Runs once per
## cohort x region, over all prepared chromosomes.
##
## Usage:
##   Rscript 02_summarize.R --cohort AA --region dlpfc --run-id ID
##
## Produces the immutable vmr_set_id that every downstream module cites.
## AGENTS.md 6: "VMR turnover never authorizes reuse of downstream numbers."

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(bsseq)
    library(HDF5Array)
    library(data.table)
    library(GenomicRanges)
})

opts <- parse_v2_args(require = c("cohort", "region", "run_id"))
cohort <- opts$cohort; region <- opts$region

th <- load_config("thresholds")
assert_locked(list(thresholds = th, cohorts = load_config("cohorts")),
              allow_unlocked = opts$allow_unlocked)
vmr_cfg <- th$vmr

## regionFinder3() is not exported by bsseq (checked against 1.42.0), so there
## is no supported call to migrate to. AGENTS.md 7.1 forbids unpinned
## `package:::` internals -- pinning here means asserting the version we
## validated against and recording it in the run manifest, so a silent upstream
## change to an internal function cannot alter VMR boundaries unnoticed.
BSSEQ_PINNED <- "1.42.0"
bsseq_version <- as.character(packageVersion("bsseq"))
if (bsseq_version != BSSEQ_PINNED) {
    stop("bsseq ", bsseq_version, " is installed but VMR calling is validated ",
         "against ", BSSEQ_PINNED, ".\n  02_summarize.R uses the unexported ",
         "bsseq:::regionFinder3(); an unverified version may change VMR ",
         "boundaries silently. Re-validate, then update BSSEQ_PINNED.")
}
if (!exists("regionFinder3", envir = asNamespace("bsseq"))) {
    stop("bsseq:::regionFinder3 not found in bsseq ", bsseq_version)
}
regionFinder3 <- get("regionFinder3", envir = asNamespace("bsseq"))

module_root <- file.path(V2_ROOT, "01_vmr_catalog")
run_dir <- file.path(module_root, "_m", "runs", opts$run_id)
vmr_dir <- file.path(run_dir, "vmr")
dir.create(vmr_dir, recursive = TRUE, showWarnings = FALSE)

chroms <- chrom_order()

## ---------------------------------------------------------------- functions

#' Call VMRs on one chromosome from residual SD.
#'
#' Runs of CpGs whose residual SD exceeds the per-chromosome quantile cutoff,
#' merged within max_gap bp, keeping runs with more than min_cpgs CpGs.
call_vmrs <- function(res_var, chrom, sd_quantile, max_gap, min_cpgs) {
    v <- res_var[order(pos)]
    sd_cut <- quantile(v$sd, probs = sd_quantile, na.rm = TRUE)
    is_high <- as.integer(!is.na(v$sd) & v$sd > sd_cut)
    if (sum(is_high) == 0) {
        return(list(vmr = data.table(chr = character(), start = integer(),
                                     end = integer(), n = integer()),
                    sd_cut = sd_cut))
    }
    found <- regionFinder3(is_high, rep(paste0("chr", chrom), nrow(v)),
                           v$pos, maxGap = max_gap, verbose = FALSE)$up
    if (is.null(found) || nrow(found) == 0) {
        return(list(vmr = data.table(chr = character(), start = integer(),
                                     end = integer(), n = integer()),
                    sd_cut = sd_cut))
    }
    found <- as.data.table(found)
    vmr <- found[n > min_cpgs, .(chr = as.character(chr),
                                 start = as.integer(start),
                                 end = as.integer(end),
                                 n = as.integer(n))]
    list(vmr = vmr, sd_cut = sd_cut)
}

## ------------------------------------------------------- part 1: call VMRs

vmr_parts <- list(); membership_parts <- list(); cutoff_parts <- list()
completed <- character(); qc_failed <- character()

for (chrom in chroms) {
    var_file <- file.path(run_dir, "pca", paste0("chr_", chrom), "res_var_all.tsv")
    if (!file.exists(var_file)) {
        qc_failed <- c(qc_failed, chrom)
        message("[vmr] chr", chrom, ": no res_var_all.tsv -- marking QC-failed")
        next
    }
    res_var <- fread(var_file)
    out <- call_vmrs(res_var, chrom, vmr_cfg$sd_quantile, vmr_cfg$max_gap,
                     vmr_cfg$min_cpgs)

    cutoff_parts[[chrom]] <- data.table(
        chr = paste0("chr", chrom), sd_cutoff = as.numeric(out$sd_cut),
        n_cpgs_tested = nrow(res_var), n_vmrs = nrow(out$vmr))

    if (nrow(out$vmr) > 0) {
        vmr_parts[[chrom]] <- out$vmr
        ## CpG -> VMR membership, needed by 05_cpg_meqtl_burden.
        gr_vmr <- GRanges(out$vmr$chr, IRanges(out$vmr$start, out$vmr$end))
        gr_cpg <- GRanges(paste0("chr", chrom),
                          IRanges(res_var$pos, width = 1))
        hits <- findOverlaps(gr_cpg, gr_vmr)
        membership_parts[[chrom]] <- data.table(
            chr = paste0("chr", chrom),
            cpg_pos = res_var$pos[queryHits(hits)],
            vmr_start = out$vmr$start[subjectHits(hits)],
            vmr_end = out$vmr$end[subjectHits(hits)],
            cpg_residual_sd = res_var$sd[queryHits(hits)])
    }
    completed <- c(completed, chrom)
    message("[vmr] chr", chrom, ": ", nrow(out$vmr), " VMRs (SD cutoff ",
            signif(out$sd_cut, 4), ")")
}

reconcile(expected = chroms, completed = completed, qc_failed = qc_failed,
          run = list(dir = run_dir), allow_failures = opts$allow_unlocked)

if (length(vmr_parts) == 0) stop("No VMRs called on any chromosome")

vmr <- sort_genomic(rbindlist(vmr_parts), chr_col = "chr", pos_col = "start")
membership <- sort_genomic(rbindlist(membership_parts), chr_col = "chr",
                           pos_col = "cpg_pos")

## Stable VMR identifiers. Derived from coordinates, so the same region in two
## runs gets the same ID and cross-run comparison is a join, not a guess.
vmr[, vmr_id := paste0(chr, ":", start, "-", end)]
membership[, vmr_id := paste0(chr, ":", vmr_start, "-", vmr_end)]

## vmr_set_id: one checksum identifying this exact catalog (AGENTS.md 9).
vmr_set_id <- paste0(
    "vmrset-", cohort, "-", region, "-",
    substr(donor_checksum(paste(vmr$vmr_id, collapse = ",")), 1, 12))
vmr[, vmr_set_id := vmr_set_id]

write_atomic(vmr[, .(chr, start, end)], file.path(vmr_dir, "vmr.bed"),
             col.names = FALSE)
write_atomic(vmr, file.path(vmr_dir, "vmr_catalog.tsv"))
write_atomic(membership, file.path(vmr_dir, "cpg_vmr_membership.tsv"))
write_atomic(rbindlist(cutoff_parts), file.path(vmr_dir, "sd_cutoffs.tsv"))

message("[vmr] ", nrow(vmr), " VMRs total | vmr_set_id = ", vmr_set_id)

## AGENTS.md 9 requires vmr_set_id in the run manifest: that is the file every
## downstream module reads to confirm which catalog it consumed. 00_new_run.R
## creates the field empty because the id does not exist until the catalog is
## called, and nothing was filling it in afterwards.
append_manifest(list(dir = run_dir),
                list(vmr_set_id = vmr_set_id, n_vmrs = nrow(vmr)))

## ------------------------------- part 2: per-VMR mean methylation phenotypes

pheno_dir <- file.path(run_dir, "vmr", "phenotypes")
dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)

arm <- cohort_def(cohort)
samples <- read_psam(arm$psam)
donor_ref <- NULL
n_written <- 0L

for (chrom in names(vmr_parts)) {
    load(file.path(run_dir, "cpg", paste0("chr_", chrom), "stats.rda"))  # BSobj
    ids <- as.character(colData(BSobj)$brnum)
    if (is.null(donor_ref)) {
        donor_ref <- ids
    } else if (!identical(donor_ref, ids)) {
        stop("Donor set differs between chromosomes: chr", chrom,
             " has ", length(ids), " donors against ", length(donor_ref),
             " elsewhere. VMR phenotypes must share one donor order.")
    }

    regions <- vmr_parts[[chrom]]
    gr <- GRanges(regions$chr, IRanges(regions$start, regions$end))
    ## One getMeth() call for the whole chromosome rather than one per VMR.
    meth_reg <- getMeth(BSobj, regions = gr, what = "perRegion")

    a <- align_by_id(data.table(FID = ids), samples,
                     id_x = ids, id_y = samples$FID, ids = ids)
    for (i in seq_len(nrow(regions))) {
        out <- data.table(FID = ids, IID = a$y$IID,
                          pheno = as.numeric(meth_reg[i, ]))
        write_atomic(out, file.path(
            pheno_dir, paste0(regions$chr[i], "_", regions$start[i], "_",
                              regions$end[i], "_meth.phen")),
            col.names = FALSE)
        n_written <- n_written + 1L
    }
    message("[pheno] chr", chrom, ": ", nrow(regions), " VMR phenotypes")
    rm(BSobj, meth_reg); gc()
}

if (n_written != nrow(vmr)) {
    stop("Wrote ", n_written, " VMR phenotype files but the catalog has ",
         nrow(vmr), " VMRs.")
}

message("[pheno] ", n_written, " VMR phenotype files written")

## PLINK --keep list for step_4.sh, in the same donor order as the phenotypes.
## Writing it here rather than in the launcher means the genotype extraction can
## only ever see the donors the catalog was actually built on.
keep <- align_by_id(data.table(FID = donor_ref), samples,
                    id_x = donor_ref, id_y = samples$FID, ids = donor_ref)
write_atomic(data.table(FID = donor_ref, IID = keep$y$IID),
             file.path(vmr_dir, "donors_plink.txt"), col.names = FALSE)

write_atomic(
    data.table(
        field = c("cohort", "region", "vmr_set_id", "n_vmrs", "n_donors",
                  "n_cpgs_in_vmrs", "donor_checksum", "chromosomes",
                  "bsseq_version", "sd_quantile", "max_gap", "min_cpgs"),
        value = c(cohort, region, vmr_set_id, nrow(vmr), length(donor_ref),
                  nrow(membership), donor_checksum(donor_ref),
                  paste(names(vmr_parts), collapse = ","),
                  bsseq_version, vmr_cfg$sd_quantile, vmr_cfg$max_gap,
                  vmr_cfg$min_cpgs)),
    file.path(vmr_dir, "summarize_summary.tsv"))

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
