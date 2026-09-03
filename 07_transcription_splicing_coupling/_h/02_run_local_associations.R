#!/usr/bin/env Rscript
#### 07_transcription_splicing_coupling -- pair-level local associations ####
##
## Usage:
##   Rscript _h/02_run_local_associations.R --run-id <id> --modality psi
##
## Fits, for every VMR-feature link in the tested universe,
##
##     feature ~ VMR_methylation + covariates
##
## and keeps the methylation coefficient's test. This recomputes the legacy
## `run_association_pass()` on the ACCEPTED Module 01 VMR boundaries; the legacy
## pair table cannot be reused, because a VMR's methylation summary is a
## function of its boundary and the boundaries were corrected.
##
## Implementation note. There are ~2.7M pairs across the three modalities, so a
## per-pair lm() call is not viable. Because every pair in a modality shares the
## same donors and the same covariate matrix, the t-statistic for the
## methylation coefficient equals the partial-correlation t-statistic after
## projecting both sides onto the covariate null space. So the covariates are
## residualised out ONCE with a QR decomposition, both sides are scaled to unit
## norm, and each pair reduces to a dot product. This is algebraically identical
## to fitting the full model per pair, not an approximation.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE", file.path(Sys.getenv("V2_REPO_ROOT", "."), "07_transcription_splicing_coupling", "_h")), "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(SummarizedExperiment)
})

MODULE <- "07_transcription_splicing_coupling"

opts <- parse_v2_args(require = c("run_id", "modality"))
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
modality <- opts$modality

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field, required = TRUE) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0 || is.na(v[1])) {
        if (required) stop("Manifest has no value for '", field, "'")
        return(NA_character_)
    }
    v[1]
}
cohort <- mf("cohort"); region <- mf("region")
vmr_run <- mf("upstream_vmr_catalog_run_id")

ts <- load_run_config("transcription_splicing", run_dir)
spec <- ts$modalities[[modality]]
if (is.null(spec)) stop("Unknown modality: ", modality)
assay_kind <- spec$assay

links_f <- file.path(run_dir, "links", paste0(modality, "-links.tsv.gz"))
if (!file.exists(links_f)) stop("No link table for ", modality, "; run stage 01")
links <- fread(links_f)
n_pairs_declared <- nrow(links)
n_features_declared <- uniqueN(links$feature_id)
message("[07] ", modality, ": ", nrow(links), " declared pairs")

## --------------------------------------------------------------- assay side
rse_f <- file.path(repo_root(), ts$assay_files[[assay_kind]][[region]])
env <- new.env(parent = emptyenv())
load(rse_f, envir = env)
objs <- mget(ls(env), envir = env)
rse <- objs[[which(vapply(objs, function(x)
    inherits(x, "SummarizedExperiment"), logical(1)))[1]]]

cd <- as.data.frame(colData(rse))
if (!"BrNum" %in% names(cd)) stop("RSE colData has no BrNum column")
cd$sample_id <- as.character(cd$BrNum)

## One library per donor. Ties are broken by RIN so the choice is deterministic
## and the discarded libraries are recorded rather than silently dropped.
rin <- if ("RIN" %in% names(cd)) as.numeric(cd$RIN) else rep(NA_real_, nrow(cd))
ord <- order(cd$sample_id, -replace(rin, is.na(rin), -Inf))
cd_ord <- cd[ord, , drop = FALSE]
keep_idx <- ord[!duplicated(cd_ord$sample_id)]
dropped_libs <- setdiff(seq_len(nrow(cd)), keep_idx)
if (length(dropped_libs)) {
    write_atomic(as.data.table(cd[dropped_libs, c("sample_id"), drop = FALSE]),
                 file.path(run_dir, "excluded",
                           paste0(modality, "-duplicate-libraries.tsv")))
}
rse <- rse[, keep_idx]
cd <- cd[keep_idx, , drop = FALSE]
colnames(rse) <- cd$sample_id

## Subset to the linked features BEFORE materialising a dense matrix. The PSI
## assay is ~691k events x 487 libraries; densifying it whole is ~2.7 GB and was
## enough to OOM-kill this stage. Only features that appear in the tested
## universe can ever be fitted, so nothing is lost by dropping the rest first.
full_lib <- if (assay_kind == "gene") colSums(assay(rse, 1L)) else NULL
linked_features <- unique(links$feature_id)
feat_rows <- rownames(rse) %in% linked_features
if (!any(feat_rows)) stop("No linked feature is present in the assay")
message("[07] assay rows: ", sum(feat_rows), "/", nrow(rse), " are linked")
rse <- rse[feat_rows, ]

mat <- as.matrix(assay(rse, 1L))
rownames(mat) <- rownames(rse)

if (assay_kind == "gene") {
    nrm <- ts$normalisation$gene
    ## Library size is a property of the whole library, so it is taken from the
    ## full assay rather than from the linked subset -- otherwise CPM would be
    ## normalised by a denominator that depends on which VMRs were linked.
    lib <- full_lib[colnames(mat)]
    if (any(!is.finite(lib)) || any(lib <= 0)) {
        stop("Zero or missing library size in the gene assay")
    }
    cpm <- t(t(mat) / lib) * 1e6
    keep <- rowMeans(cpm >= nrm$min_cpm) >= nrm$min_fraction_samples
    message("[07] gene filter: ", sum(keep), "/", nrow(mat), " genes pass")
    mat <- log2(cpm[keep, , drop = FALSE] + 1)
} else {
    nrm <- ts$normalisation$psi
    sds <- matrixStats::rowSds(mat, na.rm = TRUE)
    keep <- is.finite(sds) & sds >= nrm$min_sd
    message("[07] psi filter: ", sum(keep), "/", nrow(mat), " events pass")
    mat <- mat[keep, , drop = FALSE]
}

## Missing data. PSI events are frequently unquantified in a subset of donors
## (median 62% NA per event in caudate), and an NA anywhere breaks the shared
## residualisation. A feature is kept only if its missingness is within the
## configured bound; with max_na_fraction 0 that means quantified in every donor
## of the analysis set, so every pair is fitted on the same complete design.
##
## This is a real restriction of the tested universe, not a technicality: it
## biases the retained PSI set toward constitutively quantified events. The
## count is written to excluded/ and carried into the tested universe so the
## denominator behind any coupling claim stays auditable.
max_na <- as.numeric(nrm$max_na_fraction %||% 0)
na_frac <- rowMeans(is.na(mat))
drop_na <- na_frac > max_na
if (any(drop_na)) {
    write_atomic(data.table(feature_id = rownames(mat)[drop_na],
                            na_fraction = na_frac[drop_na],
                            max_na_fraction = max_na,
                            reason = "missingness_above_threshold"),
                 file.path(run_dir, "excluded",
                           paste0(modality, "-features-dropped-missingness.tsv")))
}
message("[07] missingness filter: ", sum(!drop_na), "/", length(drop_na),
        " features within max_na_fraction ", max_na)
mat <- mat[!drop_na, , drop = FALSE]
if (nrow(mat) == 0) stop("No feature survives the missingness filter")
n_features_after_missingness <- nrow(mat)

## ---------------------------------------------------------- covariate side
pheno_f <- file.path(repo_root(), "inputs", "phenotypes", "_m",
                     if (cohort == "AA") "phenotypes-AA.tsv" else "phenotypes-all.tsv")
pheno <- fread(pheno_f)
setnames(pheno, names(pheno), tolower(names(pheno)))
pheno <- pheno[tolower(region) == tolower(get("region"))]
pheno <- unique(pheno, by = "brnum")

props_f <- file.path(repo_root(), "inputs", "cell_proportions", "_m",
                     paste0("music-proportions-", region, ".tsv"))
props <- fread(props_f)
pn <- names(props)
id_col <- pn[tolower(pn) %in% c("sample_id", "brnum", "sample", "id")][1]
if (is.na(id_col)) stop("No sample id column in ", props_f)
if ("cell_type" %in% pn && "proportion" %in% pn) {
    props <- dcast(props, get(id_col) ~ cell_type, value.var = "proportion")
    setnames(props, "id_col", "sample_id", skip_absent = TRUE)
    setnames(props, names(props)[1], "sample_id")
} else {
    setnames(props, id_col, "sample_id")
}
props[, sample_id := as.character(sample_id)]
cell_cols <- setdiff(names(props), "sample_id")
if (length(cell_cols) == 0) stop("No cell-proportion columns in ", props_f)
## Proportions are compositional; the arcsine-square-root transform is the one
## the legacy analysis used and is kept for comparability.
for (cc in cell_cols) {
    props[[cc]] <- asin(sqrt(pmax(0, pmin(1, as.numeric(props[[cc]])))))
}

qc <- data.table(sample_id = cd$sample_id)
pick <- function(cands, label) {
    hit <- cands[cands %in% names(cd)][1]
    if (is.na(hit)) stop("RSE colData has no ", label, " column")
    hit
}
qc[, RIN := as.numeric(cd[[pick(c("RIN", "rin"), "RIN")]])]
qc[, mito_mapping_rate := as.numeric(cd[[pick(c("mito_mapping_rate", "mito_rate"),
                                              "mito_mapping_rate")]])]
qc[, percent_assigned := as.numeric(cd[[pick(c("percent_assigned", "pct_assigned"),
                                             "percent_assigned")]])]

covs <- merge(qc, pheno[, .(sample_id = as.character(brnum),
                            Age = as.numeric(agedeath),
                            Sex = as.factor(sex))], by = "sample_id")
covs <- merge(covs, props, by = "sample_id")

## ---------------------------------------------------------- methylation side
phen_dir <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs", vmr_run,
                      "vmr", "phenotypes")
wanted <- unique(links$vmr_id)
## Module 01 writes one phenotype file per VMR under the legacy underscore
## filename convention (chr1_134100_134201_meth.phen) while the canonical VMR
## identifier is chr1:134100-134201. Convert rather than renaming either.
phen_filename <- function(vid) {
    paste0(gsub("[:-]", "_", vid), "_meth.phen")
}
read_phen <- function(vid) {
    f <- file.path(phen_dir, phen_filename(vid))
    if (!file.exists(f)) return(NULL)
    d <- fread(f, header = FALSE, col.names = c("sample_id", "iid", "meth"))
    d[, .(sample_id = as.character(sample_id), meth = as.numeric(meth))]
}
message("[07] reading ", length(wanted), " VMR methylation phenotypes")
phen_list <- lapply(wanted, read_phen)
names(phen_list) <- wanted
missing_phen <- wanted[vapply(phen_list, is.null, logical(1))]
phen_list <- phen_list[!vapply(phen_list, is.null, logical(1))]
if (length(missing_phen)) {
    write_atomic(data.table(vmr_id = missing_phen, reason = "no_phen_file"),
                 file.path(run_dir, "excluded",
                           paste0(modality, "-vmrs-without-phenotype.tsv")))
}
if (length(phen_list) == 0) stop("No VMR methylation phenotypes could be read")

## ------------------------------------------------------- common donor set
donors <- Reduce(intersect, list(
    colnames(mat), covs$sample_id,
    Reduce(union, lapply(phen_list, function(d) d$sample_id))
))
covs <- covs[sample_id %in% donors]
covs <- covs[complete.cases(covs)]
donors <- intersect(donors, covs$sample_id)
if (length(donors) < ts$association$min_samples) {
    stop("Only ", length(donors), " donors with complete RNA, covariate and ",
         "methylation data; config requires ",
         ts$association$min_samples)
}
message("[07] ", length(donors), " donors in the common set")

setkey(covs, sample_id)
covs <- covs[J(donors)]
mat <- mat[, donors, drop = FALSE]

## Residual maker from the covariate design. Fitting per pair would repeat this
## decomposition ~2.7M times.
X <- model.matrix(~ ., data = as.data.frame(covs[, -"sample_id"]))
qr_x <- qr(X)
rank_x <- qr_x$rank
resid_of <- function(M) t(qr.resid(qr_x, t(M)))

## Rows scaled to unit norm so a pair's correlation is a dot product.
unitise <- function(M) {
    n <- sqrt(rowSums(M * M))
    n[!is.finite(n) | n == 0] <- NA_real_
    M / n
}

feat_res <- unitise(resid_of(mat))
ok_feat <- rownames(feat_res)[is.finite(rowSums(feat_res))]

meth_mat <- do.call(rbind, lapply(names(phen_list), function(v) {
    d <- phen_list[[v]]
    setkey(d, sample_id)
    d[J(donors)]$meth
}))
rownames(meth_mat) <- names(phen_list)
## Same guard on the methylation side: a VMR whose phenotype is missing for any
## donor in the analysis set cannot enter the shared residualisation.
meth_na <- rowSums(is.na(meth_mat)) > 0
if (any(meth_na)) {
    write_atomic(data.table(vmr_id = rownames(meth_mat)[meth_na],
                            reason = "missing_methylation_in_analysis_set"),
                 file.path(run_dir, "excluded",
                           paste0(modality, "-vmrs-with-missing-methylation.tsv")))
    message("[07] dropping ", sum(meth_na),
            " VMRs with incomplete methylation on the donor set")
    meth_mat <- meth_mat[!meth_na, , drop = FALSE]
}
if (nrow(meth_mat) == 0) stop("No VMR has complete methylation on the donor set")
meth_res <- unitise(resid_of(meth_mat))
ok_meth <- rownames(meth_res)[is.finite(rowSums(meth_res))]

links <- links[vmr_id %in% ok_meth & feature_id %in% ok_feat]
if (nrow(links) == 0) stop("No testable pairs remain after filtering")
message("[07] ", nrow(links), " testable pairs")

df_resid <- length(donors) - rank_x - 1L
if (df_resid < 1) stop("No residual degrees of freedom")

## One dot-product batch per VMR.
setorder(links, vmr_id)
split_idx <- split(seq_len(nrow(links)), links$vmr_id)
r_vec <- numeric(nrow(links))
for (v in names(split_idx)) {
    idx <- split_idx[[v]]
    r_vec[idx] <- as.numeric(feat_res[links$feature_id[idx], , drop = FALSE] %*%
                             meth_res[v, ])
}
r_vec[!is.finite(r_vec)] <- NA_real_
r_vec <- pmax(pmin(r_vec, 1 - 1e-12), -1 + 1e-12)

links[, r := r_vec]
links[, t := r * sqrt(df_resid) / sqrt(1 - r^2)]
links[, p := 2 * pt(-abs(t), df = df_resid)]
links[, fdr := p.adjust(p, method = ts$association$fdr_method)]
links[, `:=`(n = length(donors), df = df_resid, modality = modality,
             cohort = cohort, region = region, run_id = opts$run_id)]

fwrite(links, file.path(run_dir, "results",
                        paste0(modality, "-pair-associations.tsv.gz")), sep = "\t")

thr <- ts$association$fdr_threshold
vmr_summary <- links[, .(
    n_pairs_tested = .N,
    n_features_tested = uniqueN(feature_id),
    n_sig_fdr = sum(fdr < thr, na.rm = TRUE),
    any_sig_fdr = as.integer(any(fdr < thr, na.rm = TRUE)),
    min_p = min(p, na.rm = TRUE),
    min_fdr = min(fdr, na.rm = TRUE),
    max_abs_r = max(abs(r), na.rm = TRUE),
    min_distance_to_feature = min(distance, na.rm = TRUE)
), by = vmr_id]
vmr_summary[, `:=`(modality = modality, cohort = cohort, region = region,
                   run_id = opts$run_id)]
write_atomic(vmr_summary, file.path(run_dir, "results",
                                    paste0(modality, "-vmr-summary.tsv")))

## The realised tested universe, after expression/missingness filtering. Stage
## 01 records what was DECLARED; this records what was actually testable, and
## the gap between them is the part a reader needs in order to interpret a null.
realised <- data.table(
    modality = modality, cohort = cohort, region = region,
    run_id = opts$run_id,
    n_donors = length(donors),
    residual_df = df_resid,
    n_pairs_declared = n_pairs_declared,
    n_pairs_tested = nrow(links),
    n_features_declared = n_features_declared,
    n_features_after_filtering = n_features_after_missingness,
    n_features_tested = uniqueN(links$feature_id),
    n_vmrs_tested = uniqueN(links$vmr_id),
    max_na_fraction = max_na,
    n_sig_pairs = sum(links$fdr < thr, na.rm = TRUE),
    n_vmrs_coupled = sum(vmr_summary$any_sig_fdr),
    fdr_method = ts$association$fdr_method,
    fdr_threshold = thr)
write_atomic(realised, file.path(run_dir, "results",
                                 paste0(modality, "-realised-universe.tsv")))

message("[07] ", modality, ": ", sum(links$fdr < thr, na.rm = TRUE),
        " significant pairs, ", sum(vmr_summary$any_sig_fdr),
        "/", nrow(vmr_summary), " VMRs coupled")
