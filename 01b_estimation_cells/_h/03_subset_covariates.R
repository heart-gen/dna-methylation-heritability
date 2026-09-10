#### 01b / stage 03: restrict the catalog's covariates and phenotypes to the cell ####
##
## Strictly speaking this is optional for correctness: load_observed_locus()
## merges the covar, qcovar and phenotype tables against the locus BED's `fam`
## with all = FALSE, so a group-restricted fam already restricts the join. It is
## written anyway because a run must be readable on its own -- an auditor should
## not have to reason about an inner join in another module to learn which
## donors a cell was fit on -- and because the explicit row counts here are what
## make stage 04's reconciliation meaningful.
##
## Nothing is recomputed. Values are copied verbatim from the accepted Module 01
## run; only rows are dropped.
##
## Usage: Rscript 03_subset_covariates.R --run-id <estcell run id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

opts <- parse_v2_args(require = c("run_id"))

run_dir <- file.path(V2_ROOT, "01b_estimation_cells", "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(field) {
    v <- manifest$value[manifest$field == field]
    if (length(v) != 1L) stop("Run manifest lacks unique field: ", field)
    as.character(v[[1L]])
}

source_run <- mval("upstream_vmr_catalog")
covar_prefix <- mval("covar_prefix")
source_dir <- file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs", source_run)
if (!dir.exists(source_dir)) stop("Upstream Module 01 run is gone: ", source_dir)

donors <- fread(file.path(run_dir, "vmr", "donors_plink.txt"),
                header = FALSE, colClasses = "character")
setnames(donors, 1:2, c("FID", "IID"))
keep_key <- paste(donors$FID, donors$IID, sep = "::")

#' Copy one headerless FID/IID-keyed table, keeping only this cell's donors.
#'
#' Every donor of the cell must be present in the source. A cell donor that the
#' catalog has no covariates for would be dropped silently by the downstream
#' merge, changing n without a record -- the defect V1 class.
subset_by_donor <- function(src, dest, what) {
    dt <- fread(src, header = FALSE, colClasses = "character")
    if (ncol(dt) < 2) stop(what, " has fewer than 2 columns: ", src)
    key <- paste(dt[[1]], dt[[2]], sep = "::")
    if (anyDuplicated(key)) stop("Duplicate donors in ", src)
    missing <- setdiff(keep_key, key)
    if (length(missing)) {
        stop("Cell donors absent from ", what, " (", src, "): ",
             paste(head(sub("::.*$", "", missing), 10), collapse = ", "))
    }
    ## Reorder to the cell's donor order rather than subsetting in place, so
    ## every file this run writes shares one explicit order.
    out <- dt[match(keep_key, key)]
    write_atomic(out, dest, col.names = FALSE)
    nrow(out)
}

## ------------------------------------------------------------- covariates
##
## Chromosome directories are kept even though the catalog writes an identical
## covariate table into each: locus_io.R builds covs/chr_{N}/{prefix}.covar from
## the locus, and this module's job is to satisfy that contract, not to redesign
## it.
chrom_dirs <- list.dirs(file.path(source_dir, "covs"), recursive = FALSE,
                        full.names = FALSE)
chrom_dirs <- grep("^chr_", chrom_dirs, value = TRUE)
if (!length(chrom_dirs)) {
    stop("Upstream run has no covs/chr_* directories: ", source_dir)
}

n_covar <- integer(0)
for (cd in chrom_dirs) {
    dir.create(file.path(run_dir, "covs", cd), recursive = TRUE,
               showWarnings = FALSE)
    for (ext in c("covar", "qcovar")) {
        f <- paste0(covar_prefix, ".", ext)
        n <- subset_by_donor(
            file.path(source_dir, "covs", cd, f),
            file.path(run_dir, "covs", cd, f),
            paste0(ext, " for ", cd))
        n_covar <- c(n_covar, n)
    }
}
if (length(unique(n_covar)) != 1L) {
    stop("Covariate files disagree on donor count: ",
         paste(sort(unique(n_covar)), collapse = ", "))
}
message("[covs] ", length(chrom_dirs), " chromosomes x 2 files, ",
        n_covar[[1]], " donors")

## ------------------------------------------------------------- phenotypes
pheno_src <- file.path(source_dir, "vmr", "phenotypes")
pheno_dest <- file.path(run_dir, "vmr", "phenotypes")
dir.create(pheno_dest, recursive = TRUE, showWarnings = FALSE)

pheno_files <- list.files(pheno_src, pattern = "_meth\\.phen$")
if (!length(pheno_files)) stop("Upstream run has no VMR phenotypes: ", pheno_src)

vmr <- fread(file.path(run_dir, "vmr", "vmr.bed"), header = FALSE)
if (nrow(vmr) != length(pheno_files)) {
    stop("Catalog has ", nrow(vmr), " VMRs but ", length(pheno_files),
         " phenotype files")
}

n_pheno <- integer(length(pheno_files))
for (i in seq_along(pheno_files)) {
    n_pheno[[i]] <- subset_by_donor(
        file.path(pheno_src, pheno_files[[i]]),
        file.path(pheno_dest, pheno_files[[i]]),
        paste0("phenotype ", pheno_files[[i]]))
}
if (length(unique(n_pheno)) != 1L) {
    stop("VMR phenotype files disagree on donor count: ",
         paste(sort(unique(n_pheno)), collapse = ", "))
}
if (n_pheno[[1]] != nrow(donors)) {
    stop("Phenotype files carry ", n_pheno[[1]], " donors but the cell has ",
         nrow(donors))
}
message("[pheno] ", length(pheno_files), " VMR phenotypes, ", n_pheno[[1]],
        " donors")

append_manifest(list(dir = run_dir), list(
    n_covariate_chromosomes = length(chrom_dirs),
    n_phenotype_files = length(pheno_files)
))
