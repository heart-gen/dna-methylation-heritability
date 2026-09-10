#### 01b / stage 01: normalize the within-group PCA into covs/genotype_pcs.tsv ####
##
## plink2 --pca writes `#FID IID PC1 PC2 ...` (or `#IID PC1 ...` when the fileset
## carries single-column IDs). load_observed_locus() reads `FID IID snpPC1..k`.
## This stage does the rename in one place, checks the donors against the cell's
## own --keep list, and refuses a PCA that silently returned fewer components
## than the covariate model asks for.
##
## Usage: Rscript 01_write_genotype_pcs.R --run-id <estcell run id>

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

pc_names <- strsplit(mval("genotype_pcs"), ",", fixed = TRUE)[[1]]
pc_names <- trimws(pc_names)
if (!length(pc_names)) stop("Run manifest carries no genotype_pcs")

eigenvec_path <- file.path(run_dir, "covs", "pca", "group.eigenvec")
if (!file.exists(eigenvec_path)) stop("No PCA output at ", eigenvec_path)

ev <- fread(eigenvec_path, header = TRUE)
setnames(ev, sub("^#", "", names(ev)))

pc_cols <- grep("^PC[0-9]+$", names(ev), value = TRUE)
pc_cols <- pc_cols[order(as.integer(sub("^PC", "", pc_cols)))]
if (length(pc_cols) < length(pc_names)) {
    stop("PCA returned ", length(pc_cols), " components but the covariate ",
         "model requires ", length(pc_names), ". A cell this small cannot ",
         "support the configured genotype_pcs; that is a PI decision, not ",
         "something to work around here.")
}
pc_cols <- pc_cols[seq_along(pc_names)]

donors <- fread(file.path(run_dir, "vmr", "donors_plink.txt"),
                header = FALSE, colClasses = "character")
setnames(donors, 1:2, c("FID", "IID"))

## plink2 drops the FID column when the fileset has single-column IDs. Recover
## it from the cell's own donor list rather than guessing, so the join key is
## the same pair locus_io.R merges on.
if (!"FID" %in% names(ev)) {
    if (!"IID" %in% names(ev)) stop("PCA output has neither FID nor IID")
    ev[, FID := donors$FID[match(as.character(IID), donors$IID)]]
    if (anyNA(ev$FID)) stop("PCA output carries IIDs absent from the cell")
}

out <- ev[, c("FID", "IID", pc_cols), with = FALSE]
out[, (c("FID", "IID")) := lapply(.SD, as.character), .SDcols = c("FID", "IID")]
setnames(out, pc_cols, pc_names)

## The PCA must cover exactly this cell's donors. plink2 silently drops a donor
## with no genotype data, and a missing row would later be an inner-merge
## deletion inside locus_io.R -- the V1 class of defect.
assert_no_dups(out$FID, "FID in PCA output")
assert_present(out$FID, donors$FID, "donors in the PCA output")
extra <- setdiff(out$FID, donors$FID)
if (length(extra)) {
    stop("PCA output carries donors outside the cell: ",
         paste(head(extra, 10), collapse = ", "))
}
out <- out[match(donors$FID, out$FID)]
stopifnot(identical(out$FID, donors$FID))

for (nm in pc_names) {
    out[[nm]] <- as.numeric(out[[nm]])
    if (anyNA(out[[nm]])) stop("Missing values in ", nm)
    if (stats::sd(out[[nm]]) == 0) {
        stop(nm, " has zero variance within this cell; it cannot be a covariate")
    }
}

write_atomic(out, file.path(run_dir, "covs", "genotype_pcs.tsv"))
message("[pca] ", nrow(out), " donors x ", length(pc_names), " within-group PCs")
