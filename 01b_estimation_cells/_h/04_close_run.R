#### 01b / stage 04: reconcile the extraction and seal the cell ####
##
## AGENTS.md 9: "Production runs have zero tolerance for unexplained
## computational failures", and a SLURM array that exits 0 on every task is not
## evidence that every task produced output. This stage proves that every VMR in
## the catalog was accounted for in THIS cell, as one of:
##
##   extracted        a BED exists
##   no_cis_variants  a .no-snps marker exists (an explained exclusion)
##
## and nothing else. A window can carry variants in one donor group and none in
## the other, so the exclusion set is genuinely per cell and is recorded here
## rather than inherited from the pooled run.
##
## Usage: Rscript 04_close_run.R --run-id <estcell run id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

opts <- parse_v2_args(require = c("run_id"))

module_root <- file.path(V2_ROOT, "01b_estimation_cells")
run_dir <- file.path(module_root, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(field) {
    v <- manifest$value[manifest$field == field]
    if (length(v) != 1L) stop("Run manifest lacks unique field: ", field)
    as.character(v[[1L]])
}
group <- mval("estimation_group")
source_run <- mval("upstream_vmr_catalog")

## The locus set must still be the catalog's, byte for byte. If this run's
## vmr.bed ever diverged from its source, the two donor groups would no longer
## be evaluated on one fixed locus set and the whole design would be void.
source_vmr <- file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs", source_run,
                        "vmr")
cell_bed <- file.path(run_dir, "vmr", "vmr.bed")
for (f in c("vmr.bed", "vmr_catalog.tsv")) {
    if (!identical(file_sha256(file.path(source_vmr, f)),
                   file_sha256(file.path(run_dir, "vmr", f)))) {
        stop("This cell's ", f, " differs from the catalog's. The locus set ",
             "must be shared across donor groups; refusing to seal.")
    }
}

vmr <- fread(cell_bed, header = FALSE)
setnames(vmr, 1:3, c("chr", "start", "end"))
vmr[, task_id := .I]

stem <- function(i) paste0(vmr$start[i], "_", vmr$end[i])
chr_dir <- function(i) paste0("chr_", sub("^chr", "", vmr$chr[i]))
prefix <- file.path(run_dir, "plink_format", chr_dir(seq_len(nrow(vmr))),
                    paste0("TOPMed_LIBD-", group, ".",
                           stem(seq_len(nrow(vmr)))))

has_bed <- file.exists(paste0(prefix, ".bed"))
has_marker <- file.exists(paste0(prefix, ".no-snps"))

## Both present is not "extracted" -- it is a stale file from a resubmitted task
## and means the directory no longer describes one run.
both <- which(has_bed & has_marker)
if (length(both)) {
    stop("VMRs carry BOTH a BED and a .no-snps marker (stale output from a ",
         "resubmission): ",
         paste(head(vmr$task_id[both], 10), collapse = ", "))
}

completed <- vmr$task_id[has_bed]
excluded  <- vmr$task_id[has_marker]
missing   <- vmr$task_id[!has_bed & !has_marker]

## The extraction log is the array's own account of itself. Cross-check it
## against the filesystem: a task that logged "extracted" but left no BED is a
## lost file, not an exclusion.
log_dir <- file.path(run_dir, "plink_format", "extraction_log")
if (dir.exists(log_dir)) {
    log_files <- list.files(log_dir, pattern = "^task_[0-9]+\\.tsv$",
                            full.names = TRUE)
    if (length(log_files)) {
        logged <- rbindlist(lapply(log_files, fread, header = FALSE),
                            use.names = FALSE)
        setnames(logged, 1:9, c("task_id", "chr", "start", "end", "window_start",
                                "window_end", "clamped", "status", "n_variants"))
        write_atomic(logged[order(task_id)],
                     file.path(run_dir, "extraction-log.tsv"))
        claimed <- logged$task_id[logged$status == "extracted"]
        lost <- setdiff(claimed, completed)
        if (length(lost)) {
            stop("Tasks logged an extraction but left no BED: ",
                 paste(head(lost, 10), collapse = ", "))
        }
    }
}

reconcile(
    expected  = vmr$task_id,
    completed = completed,
    excluded  = excluded,
    failed    = missing,
    run = list(dir = run_dir))

if (length(missing)) {
    stop(length(missing), " VMR(s) produced neither a BED nor a .no-snps ",
         "marker. Resubmit those array tasks; do not seal a partial cell.")
}

## The genotype PCs must exist and cover the cell exactly. Sealing without them
## would let Module 02 fall back to the no-PC covariate model silently.
pc_path <- file.path(run_dir, "covs", "genotype_pcs.tsv")
if (!file.exists(pc_path)) {
    stop("No covs/genotype_pcs.tsv. Run step_1_group_pca.sh before closing; a ",
         "cell without within-group PCs is not the PI-locked covariate model ",
         "(config/covariates.yml: estimation_cells.genotype_pcs).")
}
pcs <- fread(pc_path, colClasses = list(character = c("FID", "IID")))
donors <- fread(file.path(run_dir, "vmr", "donors_plink.txt"),
                header = FALSE, colClasses = "character")
setnames(donors, 1:2, c("FID", "IID"))
if (!identical(pcs$FID, donors$FID)) {
    stop("genotype_pcs.tsv donors do not match the cell's donor list in order")
}

run <- list(run_id = opts$run_id, dir = run_dir, module_root = module_root)
append_manifest(run, list(
    n_vmrs                = nrow(vmr),
    n_vmrs_extracted      = length(completed),
    n_vmrs_no_cis_variant = length(excluded),
    n_genotype_pcs        = sum(grepl("^snpPC[0-9]+$", names(pcs))),
    catalog_vmr_bed_sha256 = file_sha256(cell_bed)
))

close_run(run)

message("[close] cell ", opts$run_id, " is sealed.\n",
        "  extracted:        ", length(completed), "\n",
        "  no cis variants:  ", length(excluded), "\n",
        "  Record it in 01b_estimation_cells/README.md with its acceptance ",
        "gate before Module 02 or 03 consumes it (AGENTS.md 6).")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
