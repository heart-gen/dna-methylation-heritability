#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- per-VMR annotation and covariate table ####
##
## Usage:
##   Rscript _h/01_build_features.R --run-id rra-AA-caudate-20260823
##
## Builds one row per interpretable VMR: the outcomes (H3K9me3 overlap,
## quiescent-chromatin overlap, LINE/L1 overlap), the primary predictor
## (standardized local SNP contribution score), the secondary predictor
## (r2_pred_oof), and every adjustment covariate named in
## config/repeat_annotations.yml.
##
## The covariate list is not decorative. VMR length, CpG count and density,
## mappability, segmental-duplication overlap and tested-SNP count all correlate
## with BOTH repeat content and the ability to estimate local genetic variance
## at all. An unadjusted overlap test recovers those confounds, not biology.
## This script therefore FAILS if a declared covariate cannot be constructed,
## rather than dropping it and proceeding with a thinner model -- and `required`
## below is generated FROM the config list, so a covariate cannot be declared in
## config and quietly omitted here.
##
## Builds: VMRs and every project asset are hg38. The Roadmap H3K9me3 and
## ChromHMM tracks are hg19, so the chromatin outcomes are computed on lifted
## coordinates. See _h/annotation_io.R.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(rtracklayer)
})

MODULE <- "04_repeat_repressive_architecture"
source(file.path(repo_root(), MODULE, "_h", "annotation_io.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mval("cohort"); region <- mval("region")
annot <- load_config("repeat_annotations")

## -------------------------------------------------------- primary predictor
lcg <- load_local_genetic_control(
    mval("upstream_local_genetic_variance_run_id"),
    region = region, cohort = cohort, eligible_only = TRUE
)

## A smoke run subsamples the catalog evenly across chromosomes. The expensive
## covariates -- WGBS coverage (one BSseq object per chromosome) and the
## cell-composition R2 (one phenotype file per VMR) -- scale with the catalog,
## so the full build is a compute-node job, not something to exercise inline.
if (!is.null(opts$smoke_n)) {
    if (!identical(mval("smoke_run"), "TRUE")) {
        stop("--smoke-n is only valid for a run opened with --allow-unlocked")
    }
    keep <- unique(round(seq(1, nrow(lcg), length.out = as.integer(opts$smoke_n))))
    lcg <- lcg[keep]
    message("[04] SMOKE: ", nrow(lcg), " VMRs sampled from the catalog")
}

vmr <- GRanges(seqnames = lcg$chrom,
               ranges = IRanges(start = lcg$start, end = lcg$end),
               vmr_id = lcg$vmr_id)

feat <- data.table(
    vmr_id = lcg$vmr_id, chrom = lcg$chrom,
    start = lcg$start, end = lcg$end,
    local_snp_contribution_score = lcg$local_snp_contribution_score,
    local_snp_contribution_score_z = lcg$local_snp_contribution_score_z,
    local_snp_contribution_quartile = lcg$local_snp_contribution_quartile
)

## ------------------------------------------------------ secondary predictor
pred_run <- mval("upstream_local_snp_prediction_run_id")
if (!is.na(pred_run) && nzchar(pred_run)) {
    pf <- list.files(file.path(repo_root(), "03_local_snp_prediction", "_m",
                               "runs", pred_run, "results", "combined"),
                     pattern = "^oof-prediction-.*\\.tsv$", full.names = TRUE)
    if (length(pf) != 1) stop("Expected one OOF prediction table for ", pred_run)
    pr <- fread(pf[1])
    if ("r_squared_cv" %in% names(pr)) {
        stop("Prediction table carries the banned legacy metric r_squared_cv ",
             "(AGENTS.md 3)")
    }
    feat <- merge(feat, pr[, .(vmr_id, r2_pred_oof)], by = "vmr_id", all.x = TRUE)
    feat[, r2_pred_oof_z := as.numeric(scale(r2_pred_oof))]
}
setkey(feat, vmr_id)
feat <- feat[lcg$vmr_id]            # restore catalog order after the merge

## ------------------------------------------------ outcomes (hg19 chromatin)
lo <- lift_vmrs_to_hg19(vmr, annot)
write_atomic(lo$report, file.path(run_dir, "results", "liftover-report.tsv"))
n_drop <- sum(lo$report$liftover_status != "unique")
message("[04] liftover hg38->hg19: ", nrow(lo$report) - n_drop, " unique, ",
        n_drop, " dropped (", sum(lo$report$liftover_status == "unmapped"),
        " unmapped, ", sum(lo$report$liftover_status == "multi_mapping"),
        " multi-mapping)")

## Every hg19 chromatin track, in one table. The first two are the primary
## repressive outcomes (BH family); the remaining four are the prespecified
## controls -- H3K27me3/bivalent for Polycomb specificity, accessible/H3K27ac
## for the complementary contrast the concentration claim needs. All six come
## from the same EIDs and the same liftover, so nothing but the annotation
## differs between them (config/repeat_annotations.yml, multiple_testing).
CHROMATIN_TRACKS <- list(
    list(key = "h3k9me3",    reader = "gappedpeak"),
    list(key = "quiescent",  reader = "chromhmm"),
    list(key = "h3k27me3",   reader = "gappedpeak"),
    list(key = "bivalent",   reader = "chromhmm"),
    list(key = "accessible", reader = "chromhmm"),
    list(key = "h3k27ac",    reader = "gappedpeak")
)

chrom_feat <- do.call(cbind, c(
    list(data.table(vmr_id = mcols(lo$gr)$vmr_id)),
    lapply(CHROMATIN_TRACKS, function(tr) {
        gr <- if (tr$reader == "gappedpeak") {
            load_gappedpeak_hg19(annot, region, tr$key)
        } else {
            load_chromhmm_states_hg19(annot, region, tr$key)
        }
        message("[04] ", tr$key, ": ", length(gr), " hg19 intervals")
        overlap_features(lo$gr, gr, tr$key)
    })
))
feat <- merge(feat, chrom_feat, by = "vmr_id", all.x = TRUE, sort = FALSE)
## A VMR that did not lift has NO chromatin outcome. It is kept in the table
## with NA and excluded by the model's complete-case rule, so the count of loci
## actually tested is visible rather than implied.
feat[, chromatin_outcome_available := !is.na(h3k9me3_frac)]

## LINE/L1 is an hg38 asset -- no liftover, and the count of tested loci differs
## from the chromatin outcomes for that reason. The gate stage reports both.
feat <- cbind(feat, overlap_features_bed(vmr, annot$repeatmasker$line_l1_bed, "line_l1"))

## ---------------------------------------------------------------- covariates
feat[, `:=`(
    vmr_length = as.numeric(end - start + 1L),
    cpg_count = as.numeric(lcg$n_cpgs),
    mean_methylation = lcg$mean_methylation,
    methylation_variance = lcg$methylation_variance,
    tested_snp_count = as.numeric(lcg$n_variants)
)]
feat[, cpg_density := cpg_count / vmr_length]

## GC content, hg38 reference sequence.
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
gc <- letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, vmr),
                      letters = "GC", as.prob = TRUE)
feat[, gc_content := as.numeric(gc)]

## Segmental duplications, problematic regions (ENCODE blacklist), and
## SNP proximity. All hg38 BED assets, all registered in the annotation manifest.
feat <- cbind(feat,
    overlap_features_bed(vmr, annot$sensitivities$exclude_segdups$bed, "segdup"),
    overlap_features_bed(vmr, annot$covariate_sources$problematic_region_overlap,
                         "problematic"),
    overlap_features_bed(vmr, annot$covariate_sources$snp_proximity, "snp_proximal"))

## Mappability: mean Umap k24 multi-read score over the VMR.
umap_p <- annot_path(annot$sensitivities$high_mappability$umap_k24_bw)
feat[, mappability := {
    sc <- summary(BigWigFile(umap_p), which = vmr, type = "mean")
    as.numeric(unlist(lapply(sc, function(x) if (length(x$score)) x$score else NA_real_)))
}]

## Broad genomic annotation, by fixed precedence (promoter > exonic > intronic).
gtf_bed <- annot_path("/projects/b1213/resources/genomes/human/gencode-v47/gtf/gencode.v47.primary_assembly.annotation.bed")
feat[, broad_genomic_annotation := broad_genomic_annotation(vmr, gtf_bed)]

## Cell-composition-associated methylation. AGENTS.md 7.4 asks for
## "cell-composition-associated methylation properties" as an adjustment. The
## donor deconvolution estimates are per DONOR, not per VMR, so the per-VMR
## quantity is how strongly a VMR's methylation tracks composition: the R^2 of
## its methylation across donors on the donor cell-proportion principal
## components. A VMR whose signal is mostly composition is exactly the one whose
## repeat enrichment would otherwise be over-read.
prop_f <- file.path(repo_root(), "inputs", "cell_proportions", "_m",
                    sprintf("dnam-scmd-proportions-%s.tsv", region))
if (!file.exists(prop_f)) stop("Cell-proportion estimates not found: ", prop_f)
props <- fread(prop_f)
wide <- dcast(props, sample_id ~ cell_type, value.var = "proportion")
pmat <- as.matrix(wide[, -1]); rownames(pmat) <- wide$sample_id
pmat <- pmat[, apply(pmat, 2, function(z) is.finite(stats::var(z)) &&
                                          stats::var(z) > 0), drop = FALSE]
n_pc <- min(as.integer(annot$covariate_sources$n_cell_composition_pcs %||% 3L),
            ncol(pmat) - 1L)
pcs <- stats::prcomp(pmat, center = TRUE, scale. = TRUE)$x[, seq_len(n_pc), drop = FALSE]

vmr_run_dir <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs",
                         mval("upstream_vmr_catalog_run_id"))
phen_dir <- file.path(vmr_run_dir, "vmr", "phenotypes")
if (!dir.exists(phen_dir)) stop("Module 01 phenotypes not found: ", phen_dir)

feat[, cell_composition_r2 := {
    vapply(seq_len(.N), function(i) {
        f <- file.path(phen_dir, sprintf("%s_%d_%d_meth.phen",
                                         chrom[i], start[i], end[i]))
        if (!file.exists(f)) return(NA_real_)
        ph <- fread(f, header = FALSE, col.names = c("FID", "IID", "y"),
                    colClasses = list(character = 1:2))
        idx <- match(ph$FID, rownames(pcs))
        ok <- !is.na(idx) & is.finite(ph$y)
        if (sum(ok) < 20) return(NA_real_)
        summary(stats::lm(ph$y[ok] ~ pcs[idx[ok], , drop = FALSE]))$r.squared
    }, numeric(1))
}]

## WGBS coverage: mean per-CpG read depth over the VMR's constituent CpGs,
## read from the same BSseq objects Module 01 built the catalog from. Done once
## per chromosome, because loading a BSobj per VMR would be 11k loads.
membership <- fread(file.path(vmr_run_dir, "vmr", "cpg_vmr_membership.tsv"))
cov_by_vmr <- rbindlist(lapply(sort(unique(feat$chrom)), function(cc) {
    chrom_n <- sub("^chr", "", cc)
    bs <- tryCatch(load_bsobj(region, chrom_n, require_smoothed = FALSE),
                   error = function(e) NULL)
    if (is.null(bs)) return(NULL)
    obj <- if (is.list(bs) && !is.null(bs$bsobj)) bs$bsobj else bs
    pos <- BiocGenerics::start(SummarizedExperiment::rowRanges(obj))
    depth <- DelayedMatrixStats::rowMeans2(bsseq::getCoverage(obj, type = "Cov"))
    mem <- membership[chr == cc]
    mem[, depth := depth[match(cpg_pos, pos)]]
    mem[, .(wgbs_coverage = mean(depth, na.rm = TRUE)), by = vmr_id]
}), fill = TRUE)
feat <- merge(feat, cov_by_vmr, by = "vmr_id", all.x = TRUE, sort = FALSE)

## ------------------------------------------------------------- completeness
## The declared covariate list drives this check, so config and code cannot
## drift apart. Each declared name maps to the column(s) that realize it.
realized <- list(
    vmr_length = "vmr_length", cpg_count = "cpg_count",
    cpg_density = "cpg_density", gc_content = "gc_content",
    mean_methylation = "mean_methylation",
    methylation_variance = "methylation_variance",
    wgbs_coverage = "wgbs_coverage", tested_snp_count = "tested_snp_count",
    snp_proximity = "snp_proximal_frac", mappability = "mappability",
    segdup_overlap = "segdup_frac",
    problematic_region_overlap = "problematic_frac",
    broad_genomic_annotation = "broad_genomic_annotation",
    cell_composition_pcs = "cell_composition_r2"
)
## Both lists are built here. `descriptive_covariates` are the ones the
## 2026-09-02 amendment removed from the primary FORMULA as mediators or as
## constitutive of an outcome -- they are still constructed, still checked, and
## still fail loudly, because the sensitivity arms and the post-hoc QC ladder
## need the columns. Narrowing the model must not narrow the feature table.
declared <- unique(c(unlist(annot$covariates),
                     unlist(annot$descriptive_covariates)))
undeclared <- setdiff(declared, names(realized))
if (length(undeclared)) {
    stop("config declares covariate(s) this script does not build: ",
         paste(undeclared, collapse = ", "),
         "\n  Implement them here; do not fit a thinner model.")
}
required <- unlist(realized[declared])
missing_cols <- setdiff(required, names(feat))
if (length(missing_cols) > 0) {
    stop("Declared covariates could not be constructed: ",
         paste(missing_cols, collapse = ", "))
}
incomplete <- required[vapply(required, function(k) all(is.na(feat[[k]])), logical(1))]
if (length(incomplete) > 0) {
    stop("Covariate(s) constructed but entirely NA: ",
         paste(incomplete, collapse = ", "))
}

## Per-covariate completeness is REPORTED, not silently absorbed by the model's
## complete-case rule.
completeness <- data.table(
    covariate = required,
    n_missing = vapply(required, function(k) sum(is.na(feat[[k]])), integer(1)),
    frac_missing = vapply(required, function(k) mean(is.na(feat[[k]])), numeric(1))
)
write_atomic(completeness, file.path(run_dir, "results", "covariate-completeness.tsv"))
print(completeness)

feat[, `:=`(region = region, population = cohort, vmr_set_id = mval("vmr_set_id"))]
write_atomic(feat, file.path(run_dir, "results", "vmr-features.tsv"))
message("[04] built features for ", nrow(feat), " interpretable VMRs")
