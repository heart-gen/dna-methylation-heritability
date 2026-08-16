#### 01_vmr_catalog / 01_analyze: methylation PCs and residual variance ####
##
## Merges the legacy 02.pca.R and 02b.res_var.R into one script, because the
## split is what allowed the two steps to disagree about donor order.
##
## Usage:
##   Rscript 01_analyze.R --cohort AA --region dlpfc --chrom 21 --run-id ID
##
## THIS SCRIPT CARRIES THE V1 FIX.
##
## Legacy 02b.res_var.R:39-40 did:
##     meth_levels <- meth_levels[match(valid_ids, brain_id), , drop = FALSE]
##     pc_filt     <- pc[V1 %in% valid_ids, -1, with = FALSE]
## The response was reordered to valid_ids order; the design kept pc.csv order.
## Measured: caudate chr1 153/153 rows wrong, DLPFC 92/96, hippocampus 98/101.
## Consequence: the 99th-percentile SD cutoff shifted +9.2% and ~10% of seed
## CpGs changed, invalidating every VMR set in the repository.
##
## In v2 every pairing of a response with a design goes through align_by_id(),
## which reorders BOTH sides and asserts the IDs match. See
## 00_shared/identity.R and 00_shared/tests/test-identity.R.
##
## Other repairs:
##   V6  top-CpG indexing guarded by the observed CpG count.
##   V7  chunks recombined in numeric order with a donor-order check.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(bsseq)
    library(HDF5Array)
    library(DelayedMatrixStats)
    library(data.table)
    library(matrixStats)
    library(limma)
})

opts <- parse_v2_args(require = c("cohort", "region", "chrom", "run_id"))
cohort <- opts$cohort; region <- opts$region; chrom <- as.character(opts$chrom)

th   <- load_config("thresholds")
covs <- load_config("covariates")
assert_locked(list(thresholds = th, cohorts = load_config("cohorts")),
              allow_unlocked = opts$allow_unlocked)

arm <- cohort_def(cohort)
vmr_cfg <- th$vmr
n_snp_pcs  <- vmr_cfg$snp_pcs_regressed
n_meth_pcs <- vmr_cfg$meth_pcs_regressed

module_root <- file.path(V2_ROOT, "01_vmr_catalog")
run_dir <- file.path(module_root, "_m", "runs", opts$run_id)
base <- if (is_primary_chrom(chrom)) run_dir else file.path(run_dir, "excluded")
cpg_dir <- file.path(base, "cpg", paste0("chr_", chrom))
pca_dir <- file.path(base, "pca", paste0("chr_", chrom))
dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------- functions

#' Genotype PCs for the analysis donors, aligned to a given donor order.
#'
#' Donors with any NA among the requested PCs are dropped, loudly.
get_snp_pcs <- function(pheno_file, ids, n_pcs) {
    pc_cols <- paste0("snpPC", seq_len(n_pcs))
    tab <- fread(pheno_file, header = TRUE)
    missing_cols <- setdiff(pc_cols, names(tab))
    if (length(missing_cols) > 0) {
        stop("Phenotype table lacks genotype PC column(s): ",
             paste(missing_cols, collapse = ", "))
    }
    ## The phenotype table has one row per donor x region; collapse to donors,
    ## since genotype PCs do not vary by brain region.
    tab <- unique(tab[, c("brnum", pc_cols), with = FALSE], by = "brnum")
    tab <- tab[brnum %in% ids]

    na_donors <- tab$brnum[rowSums(is.na(tab[, pc_cols, with = FALSE])) > 0]
    if (length(na_donors) > 0) {
        message("[snpPC] dropping ", length(na_donors),
                " donor(s) with NA genotype PCs: ",
                paste(head(na_donors, 5), collapse = ", "))
        tab <- tab[!brnum %in% na_donors]
    }
    tab
}

#' Regress genotype PCs out of methylation, keeping population structure from
#' driving the methylation PCs.
residualize_on_snp_pcs <- function(meth, meth_ids, snp_pcs, n_pcs) {
    ## meth is CpG x donor; align donors (columns) to the genotype PC rows.
    a <- align_by_id(t(meth), snp_pcs, id_x = meth_ids, id_y = snp_pcs$brnum)
    design <- cbind(1, as.matrix(a$y[, paste0("snpPC", seq_len(n_pcs)),
                                     with = FALSE]))
    stopifnot(nrow(design) == nrow(a$x))
    fit <- limma::lmFit(t(a$x), design)
    list(resid = limma::residuals.MArrayLM(fit, t(a$x)), ids = a$ids)
}

#' Regress the top methylation PCs out of a CpG chunk.
#'
#' The residual SD of these values is what defines a VMR, so a misalignment
#' here propagates into the region definitions themselves.
residualize_on_meth_pcs <- function(meth_chunk, chunk_ids, pc_mat, pc_ids, n_pcs) {
    a <- align_by_id(meth_chunk, pc_mat, id_x = chunk_ids, id_y = pc_ids)
    design <- cbind(1, as.matrix(a$y[, seq_len(n_pcs), drop = FALSE]))
    stopifnot(nrow(design) == nrow(a$x))
    fit <- limma::lmFit(t(a$x), design)
    list(resid = t(limma::residuals.MArrayLM(fit, t(a$x))), ids = a$ids)
}

## ------------------------------------------------- part 1: methylation PCs

message("[load] ", file.path(cpg_dir, "stats.rda"))
load(file.path(cpg_dir, "stats.rda"))   # sds, means, BSobj

meth_ids <- as.character(colData(BSobj)$brnum)

## V6: cap by the observed CpG count. The legacy DLPFC/hippocampus copies did
## v_top[1:10^6, ], which fabricates NA rows on any chromosome with fewer than
## a million surviving CpGs.
v <- data.table(chr = chrom, start = start(BSobj), sd = sds)
v_top <- v[order(-sd)][seq_len(min(vmr_cfg$top_variable_cpgs, nrow(v)))]
message("[pca] using ", nrow(v_top), " of ", nrow(v), " CpGs (cap ",
        vmr_cfg$top_variable_cpgs, ")")

BS_top <- BSobj[is.element(start(BSobj), v_top$start), ]
meth_top <- as.matrix(getMeth(BS_top))
colnames(meth_top) <- as.character(colData(BS_top)$brnum)

snp_pcs <- get_snp_pcs(arm$phenotype_table, meth_ids, n_snp_pcs)
res <- residualize_on_snp_pcs(meth_top, colnames(meth_top), snp_pcs, n_snp_pcs)

## prcomp() on donors x CpGs. Row order of pc$x is res$ids, and we carry that
## vector forward explicitly rather than re-deriving it later -- re-deriving is
## exactly how V1 happened.
pc <- prcomp(t(res$resid), scale. = TRUE, center = TRUE)
stopifnot(identical(rownames(pc$x), res$ids))
message("[pca] cumulative variance explained (PC1-10): ",
        paste(round(cumsum(summary(pc)$importance[2, seq_len(min(10, ncol(pc$x)))]), 4),
              collapse = " "))

pc_dt <- data.table(brnum = res$ids, pc$x)
write_atomic(pc_dt, file.path(pca_dir, "pc.csv"), sep = ",")

pdf(file.path(pca_dir, "pca.pdf"))
plot(pc, main = paste0(cohort, " ", region, " chr", chrom))
plot(pc$x[, 1], pc$x[, 2], xlab = "PC1", ylab = "PC2",
     main = paste0(cohort, " ", region, " chr", chrom))
dev.off()

## --------------------------------- part 2: residual variance, chunk by chunk

meth_file <- file.path(cpg_dir, "cpg_meth.phen")
meth_rds <- file.path(cpg_dir, "cpg_meth.rds")
## Prefer the binary twin: cpg_meth.phen is donors x CpGs, and fread() segfaults
## nondeterministically on the >1M-column chromosomes (see 00_prepare.R).
use_rds <- file.exists(meth_rds)
if (use_rds) {
    meth_bin <- readRDS(meth_rds)
    cpg_cols <- meth_bin$cpg
} else {
    ## nrows = 0 makes fread fall back to V1..Vn instead of reading the header, so
    ## the column names silently become positional and every select() misses.
    header <- names(fread(meth_file, nrows = 1L, header = TRUE))
    if (!all(c("FID", "IID") %in% header)) {
        stop("cpg_meth.phen is missing FID/IID columns: ", meth_file)
    }
    cpg_cols <- setdiff(header, c("FID", "IID"))
}
n_cpg <- length(cpg_cols)
chunk_size <- vmr_cfg$chunk_cols
starts <- seq(1, n_cpg, by = chunk_size)
message("[resid] ", n_cpg, " CpGs in ", length(starts), " chunk(s) of ",
        chunk_size)

pc_mat <- as.matrix(pc_dt[, paste0("PC", seq_len(n_meth_pcs)), with = FALSE])
pc_ids <- pc_dt$brnum

res_var_parts <- vector("list", length(starts))
resid_parts   <- vector("list", length(starts))
donor_ids_ref <- NULL

## Read the matrix ONCE. It is wide but short -- ~490k CpGs x ~150 donors, a few
## hundred MB -- so one read plus in-memory chunking is far cheaper than 98
## fread(select=) passes over a 490k-column file.
##
## header = TRUE is not optional here: the ID columns are character in both the
## header row and the data rows, so fread's auto-detection can conclude there is
## no header, rename every column to V1..Vn, and silently drop a select() as
## "not found" -- leaving a chunk with no donors.
if (use_rds) {
    message("[resid] reading methylation matrix (cpg_meth.rds)")
    meth_ids_file <- as.character(meth_bin$FID)
    meth_mat <- meth_bin$meth
    colnames(meth_mat) <- meth_bin$cpg
    rm(meth_bin); gc()
} else {
    message("[resid] reading methylation matrix (cpg_meth.phen)")
    meth_all <- fread(meth_file, header = TRUE)
    if (!all(c("FID", "IID") %in% names(meth_all)) || nrow(meth_all) == 0) {
        stop("cpg_meth.phen did not return FID/IID and rows: ", meth_file)
    }
    meth_ids_file <- as.character(meth_all$FID)
    meth_mat <- as.matrix(meth_all[, -c("FID", "IID"), with = FALSE])
    rm(meth_all); gc()
}
if (nrow(meth_mat) == 0) stop("Empty methylation matrix: ", cpg_dir)
stopifnot(identical(colnames(meth_mat), cpg_cols))

for (i in seq_along(starts)) {
    from <- starts[[i]]
    to   <- min(from + chunk_size - 1L, n_cpg)

    meth_chunk <- meth_mat[, from:to, drop = FALSE]
    r <- residualize_on_meth_pcs(meth_chunk, meth_ids_file,
                                 pc_mat, pc_ids, n_meth_pcs)

    if (is.null(donor_ids_ref)) {
        donor_ids_ref <- r$ids
    } else if (!identical(donor_ids_ref, r$ids)) {
        stop("Donor set changed between chunks ", i - 1, " and ", i,
             ". Every chunk must resolve to the same donors in the same order.")
    }

    res_var_parts[[i]] <- data.table(
        chr = chrom,
        pos = as.integer(cpg_cols[from:to]),
        sd  = colSds(r$resid),
        var = colVars(r$resid))

    rp <- data.table(r$resid)
    setnames(rp, cpg_cols[from:to])
    resid_parts[[i]] <- rp

    if (i %% 10 == 0 || i == length(starts)) {
        message("[resid] chunk ", i, "/", length(starts), ": CpGs ", from, "-", to)
    }
    rm(meth_chunk, r)
}
rm(meth_mat); gc()

## V7: parts are assembled in the order they were generated -- which is numeric
## by construction here -- and the donor vector was checked identical at every
## step above, so the cbind is safe.
res_var <- rbindlist(res_var_parts)
res_var <- res_var[order(pos)]
write_atomic(res_var, file.path(pca_dir, "res_var_all.tsv"))

samples <- read_psam(arm$psam)
a <- align_by_id(data.table(FID = donor_ids_ref), samples,
                 id_x = donor_ids_ref, id_y = samples$FID, ids = donor_ids_ref)
res_meth <- cbind(data.table(FID = donor_ids_ref, IID = a$y$IID),
                  do.call(cbind, resid_parts))
write_atomic(res_meth, file.path(cpg_dir, "res_cpg_meth.phen"))

write_atomic(
    data.table(
        field = c("cohort", "region", "chrom", "n_donors", "n_cpgs",
                  "n_chunks", "donor_checksum", "snp_pcs_regressed",
                  "meth_pcs_regressed", "sd_median", "sd_q99"),
        value = c(cohort, region, chrom, length(donor_ids_ref), nrow(res_var),
                  length(starts), donor_checksum(donor_ids_ref),
                  n_snp_pcs, n_meth_pcs,
                  format(median(res_var$sd, na.rm = TRUE), digits = 10),
                  format(quantile(res_var$sd, vmr_cfg$sd_quantile, na.rm = TRUE),
                         digits = 10))),
    file.path(pca_dir, "analyze_summary.tsv"))

message("[done] chr", chrom, ": residual variance for ", nrow(res_var),
        " CpGs x ", length(donor_ids_ref), " donors")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
