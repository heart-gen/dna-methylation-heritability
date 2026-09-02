#!/usr/bin/env Rscript
##
## 06_qc_predictability_artifact.R -- is the upper mode of r2_pred_oof real?
##
## The accepted caudate run has median r2_pred_oof ~ 0 but 20.9% of loci above
## 0.5 and 3.1% above 0.9 (max 0.970). An out-of-fold R2 of 0.93 for regional
## methylation from cis SNPs in ~153 donors is not a biologically ordinary
## result, and the leakage tripwire in 03_combine_oof.R gates on the MEDIAN, so
## it is structurally blind to leakage confined to a subpopulation.
##
## Hypothesis under test: the high mode sits in segmental-duplication / CNV /
## low-mappability territory, where WGBS "methylation" partly tracks read depth
## and copy number, and copy number is near-deterministically tagged by the
## local haplotype. Long + SNP-dense + perfectly predictable is that signature.
##
## CRITICAL CONFOUND: high-r2 loci are already known to be ~5x LONGER than the
## rest (2572 bp vs 513 bp). A longer interval overlaps any annotation more
## often by chance alone, so a raw overlap comparison is uninterpretable. Every
## test here is therefore length-adjusted, and the length-matched contrast is
## the one to read.
##
## Reads the sealed run READ-ONLY (AGENTS.md 5.2: runs are immutable) and writes
## to _m/qc/{run_id}/, outside runs/.
##
## Three checks, all reproducible from sealed runs:
##   A. ANNOTATION      does the high mode sit in segdup / blacklist / low-map
##                      territory? (length-matched, because high-r2 loci are ~5x
##                      longer and would overlap anything more often by chance)
##   B. RELATEDNESS     the remaining leakage route. Duplicate donors are already
##                      excluded by construction (donor IDs are unique), but two
##                      relatives split across outer folds would inflate OOF r2 at
##                      exactly the loci where local haplotype drives methylation.
##                      KING kinship via plink2; flags any related pair that CV
##                      splits.
##   C. CONCORDANCE     genotype is IDENTICAL across brain regions for a shared
##                      donor, so a locus predictable because of a tagged variant
##                      must be predictable in the other regions too. A region-
##                      specific high tail would mean the cause is not genetic.
##                      This is what showed the caudate ">0.9 excess" to be a
##                      threshold effect of an arbitrary cutoff on a uniformly
##                      sample-size-shifted distribution, not an anomaly.
##
## Usage:
##   Rscript _h/06_qc_predictability_artifact.R --run-id lsp-AA-caudate-20260825
##   Rscript _h/06_qc_predictability_artifact.R --run-id lsp-AA-caudate-20260825 \
##       --compare-runs lsp-AA-dlpfc-20260825,lsp-AA-hippocampus-20260825
##   ... --skip-relatedness   to omit check B (plink2 not on PATH)
##
suppressPackageStartupMessages({
    library(data.table); library(GenomicRanges); library(rtracklayer)
})

MODULE <- "03_local_snp_prediction"
source_shared <- function(f) {
    root <- Sys.getenv("V2_REPO_ROOT", unset = NA)
    if (is.na(root)) {
        d <- normalizePath(getwd())
        while (d != dirname(d) && !dir.exists(file.path(d, ".git"))) d <- dirname(d)
        root <- d
    }
    source(file.path(root, "00_shared", f)); root
}
ROOT <- source_shared("config.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, default = NULL) {
    i <- match(k, args)
    if (is.na(i) || i == length(args)) return(default)
    args[i + 1L]
}
RUN_ID <- get_arg("--run-id")
if (is.null(RUN_ID)) stop("--run-id is required")
HI <- as.numeric(get_arg("--high-threshold", "0.9"))
COMPARE <- get_arg("--compare-runs")
COMPARE <- if (is.null(COMPARE)) character(0) else strsplit(COMPARE, ",")[[1]]
SKIP_REL <- "--skip-relatedness" %in% args
PLINK2 <- Sys.getenv("PLINK2", unset = "/projects/p32505/opt/bin/plink2")
## Second-degree. KING kinship: ~0.25 MZ/dup, 0.125 first-degree, 0.0625 second.
KIN_CUT <- as.numeric(get_arg("--kinship-cutoff", "0.0884"))

read_vmrs <- function(rid) {
    d <- file.path(ROOT, MODULE, "_m", "runs", rid, "results", "combined")
    f <- list.files(d, pattern = "^oof-prediction-.*-vmrs\\.tsv$", full.names = TRUE)
    if (length(f) != 1L) stop("Expected one vmrs table for ", rid, "; got ", length(f))
    fread(f[1])
}

run_dir <- file.path(ROOT, MODULE, "_m", "runs", RUN_ID)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
comb <- list.files(file.path(run_dir, "results", "combined"),
                   pattern = "^oof-prediction-.*-vmrs\\.tsv$", full.names = TRUE)
if (length(comb) != 1L) stop("Expected exactly one oof-prediction vmrs table; got ", length(comb))

dt <- fread(comb[1])
message("Loaded ", nrow(dt), " loci from ", basename(comb[1]))

## ------------------------------------------------------------------ features
dt[, length_bp := end - start]
dt[, group := fifelse(r2_pred_oof > HI, "high", "rest")]

vmr <- GRanges(dt$chrom, IRanges(dt$start + 1L, dt$end))

sf <- file.path(ROOT, "inputs", "supportfiles", "_m")
load_bed <- function(f) {
    g <- fread(cmd = paste("zcat", shQuote(file.path(sf, f))),
               header = FALSE, select = 1:3,
               col.names = c("chrom", "start", "end"))
    g <- g[grepl("^chr[0-9]+$", chrom)]
    reduce(GRanges(g$chrom, IRanges(g$start + 1L, g$end)))
}

## Fraction of each VMR covered, not a yes/no flag: a 1 bp clip of a 2.5 kb VMR
## is not the same event as full containment, and the yes/no form is what makes
## the length confound bite hardest.
covered_fraction <- function(query, subject) {
    hits <- as.data.table(findOverlaps(query, subject))
    if (nrow(hits) == 0L) return(rep(0, length(query)))
    ov <- pintersect(query[hits$queryHits], subject[hits$subjectHits])
    agg <- hits[, .(bp = sum(width(ov))), by = queryHits]
    out <- rep(0, length(query))
    out[agg$queryHits] <- pmin(1, agg$bp / width(query)[agg$queryHits])
    out
}

for (nm in c(segdup = "genomicSuperDups.hg38.bed.gz",
             blacklist = "hg38-blacklist.v2.bed.gz",
             line_l1 = "repeat-masker.LINE_L1.hg38.bed.gz")) {
    key <- names(which(c(segdup = "genomicSuperDups.hg38.bed.gz",
                         blacklist = "hg38-blacklist.v2.bed.gz",
                         line_l1 = "repeat-masker.LINE_L1.hg38.bed.gz") == nm))
    message("  overlapping ", key, " ...")
    dt[[paste0(key, "_frac")]] <- covered_fraction(vmr, load_bed(nm))
}

## Mappability: 2 GB bigwig, so read only the VMR intervals.
bw <- file.path(sf, "k24.Umap.MultiTrackMappability.bw")
message("  summarizing mappability over ", length(vmr), " intervals ...")
## summary() returns a GRangesList; unlist() gives a GRanges, so the score has
## to be pulled from mcols -- as.numeric() on the GRanges silently yields NaN.
## Ranges are also clamped to the bigwig seqlengths first, since a VMR running
## past a contig end makes the whole call warn and return nothing useful.
dt[, mappability := tryCatch({
    bwf <- BigWigFile(bw)
    si  <- seqinfo(bwf)
    q   <- vmr[as.character(seqnames(vmr)) %in% seqlevels(si)]
    seqlevels(q) <- seqlevels(si); seqinfo(q) <- si
    q   <- trim(q)
    sc  <- mcols(unlist(summary(bwf, q, type = "mean")))$score
    out <- rep(NA_real_, length(vmr))
    out[as.character(seqnames(vmr)) %in% seqlevels(si)] <- as.numeric(sc)
    out
}, error = function(e) { warning("mappability failed: ", conditionMessage(e)); NA_real_ })]

## ------------------------------------------------------------------- testing
## Length-adjusted. A raw overlap test would be confounded by the 5x length
## difference between the groups and is reported only for contrast.
dt[, is_high := as.integer(group == "high")]
dt[, log_len := log10(length_bp)]

fit_one <- function(feature) {
    f <- dt[[feature]]
    if (all(is.na(f))) return(NULL)
    raw <- tryCatch(coef(summary(glm(is_high ~ f, family = binomial, data = dt)))["f", ],
                    error = function(e) rep(NA_real_, 4))
    adj <- tryCatch(coef(summary(glm(is_high ~ f + log_len + median_n_variants,
                                     family = binomial, data = dt)))["f", ],
                    error = function(e) rep(NA_real_, 4))
    data.table(feature = feature,
               mean_high = mean(f[dt$is_high == 1], na.rm = TRUE),
               mean_rest = mean(f[dt$is_high == 0], na.rm = TRUE),
               or_raw = exp(raw[1]), p_raw = raw[4],
               or_adj_len_variants = exp(adj[1]), p_adj = adj[4])
}
tests <- rbindlist(lapply(c("segdup_frac", "blacklist_frac", "line_l1_frac",
                            "mappability"), fit_one), fill = TRUE)

## Length-matched contrast: for each high locus take the nearest-length control
## from `rest`, sampled without replacement. This is the headline comparison.
set.seed(20260902)
hi_idx <- which(dt$is_high == 1); lo_idx <- which(dt$is_high == 0)
pool <- data.table(idx = lo_idx, len = dt$length_bp[lo_idx], used = FALSE)
setkey(pool, len)
matched <- integer(0)
for (i in hi_idx) {
    avail <- pool[used == FALSE]
    if (nrow(avail) == 0L) break
    j <- avail[which.min(abs(len - dt$length_bp[i])), idx]
    matched <- c(matched, j); pool[idx == j, used := TRUE]
}
mt <- rbindlist(lapply(c("segdup_frac", "blacklist_frac", "line_l1_frac",
                         "mappability"), function(f) {
    a <- dt[[f]][hi_idx]; b <- dt[[f]][matched]
    data.table(feature = f, mean_high = mean(a, na.rm = TRUE),
               mean_matched = mean(b, na.rm = TRUE),
               p_wilcox = tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA_real_))
}))

bins <- dt[, .(n = .N, mean_len = mean(length_bp),
               segdup = mean(segdup_frac), blacklist = mean(blacklist_frac),
               line_l1 = mean(line_l1_frac), mappability = mean(mappability, na.rm = TRUE)),
           by = .(bin = cut(r2_pred_oof, c(-Inf, .05, .5, .9, Inf),
                            labels = c("<=0.05", "0.05-0.5", "0.5-0.9", ">0.9")))][order(bin)]

## ------------------------------------------------------- B. relatedness (CV)
## Duplicate donors are already impossible here (donor IDs are unique per run),
## so the live risk is CRYPTIC relatedness: two relatives under different brain
## IDs, split across outer folds. That leaks haplotype, and it leaks hardest at
## precisely the loci this script is about. Reported as a fold-splitting count,
## not just a pair count -- a related pair kept together in every fold does not
## leak, and a pair split in one repeat out of five leaks only in that repeat.
rel_pairs <- data.table(); rel_summary <- data.table(); rel_max <- NA_real_
if (!SKIP_REL) {
    folds  <- fread(file.path(run_dir, "donor-folds.tsv"))
    donors <- unique(folds$donor)
    tmp    <- file.path(tempdir(), paste0("kin_", RUN_ID))
    keep_f <- paste0(tmp, ".keep")
    fwrite(data.table(FID = sub("::.*", "", donors), IID = sub(".*::", "", donors)),
           keep_f, sep = "\t", col.names = FALSE)

    gt <- file.path(ROOT, "inputs", "genotypes", "TOPMed_LIBD.AA")
    if (!file.exists(paste0(gt, ".pgen"))) stop("Genotypes not found: ", gt, ".pgen")
    if (!file.exists(PLINK2)) stop("plink2 not found at ", PLINK2,
                                   " (set PLINK2=, or pass --skip-relatedness)")

    ## The shipped .psam is HEADERLESS, which plink2 rejects with
    ##   "Line 1 of ....psam has fewer tokens than expected".
    ## Repair it into a temp copy rather than touching the shared input, and
    ## address the triple explicitly (--pgen/--pvar/--psam) so the repaired psam
    ## is the one used. Without this the calls below fail and -- before this was
    ## fixed -- the empty result was reported as "no related pairs", turning a
    ## broken check into a reassuring one.
    psam_fixed <- paste0(tmp, ".psam")
    src <- fread(paste0(gt, ".psam"), header = FALSE)
    if (startsWith(as.character(src[[1]][1]), "#")) {
        file.copy(paste0(gt, ".psam"), psam_fixed, overwrite = TRUE)
    } else {
        out <- src[, 1:3]; setnames(out, c("#FID", "IID", "SEX"))
        ## na = "NA" is load-bearing: fwrite's default na = "" writes the SEX
        ## column as an empty field, and plink2 then rejects the file with
        ## "Line 2 ... has fewer tokens than expected" -- the same shape of
        ## error as the missing header this block exists to repair.
        fwrite(out, psam_fixed, sep = "\t", na = "NA", quote = FALSE)
    }
    trio <- c("--pgen", paste0(gt, ".pgen"), "--pvar", paste0(gt, ".pvar"),
              "--psam", psam_fixed)

    run_plink <- function(a, what) {
        log <- paste0(tmp, ".", what, ".log")
        st  <- system2(PLINK2, a, stdout = log, stderr = log)
        if (!identical(as.integer(st), 0L))
            stop("plink2 ", what, " failed (exit ", st, "). Log:\n",
                 paste(readLines(log, warn = FALSE), collapse = "\n"))
        invisible(st)
    }
    ## KING on correlated markers is biased, so LD-prune first.
    run_plink(c(trio, "--keep", keep_f, "--autosome", "--maf", "0.05",
                "--geno", "0.05", "--indep-pairwise", "200", "50", "0.1",
                "--out", tmp), "prune")
    if (!file.exists(paste0(tmp, ".prune.in"))) stop("plink2 produced no prune.in")

    ## Filter at -1 so EVERY pair is returned: we want the observed maximum on
    ## record, not just the pairs above the cutoff. "Max kinship 0.04" is
    ## evidence; an empty table is indistinguishable from a failed run.
    run_plink(c(trio, "--keep", keep_f, "--autosome",
                "--extract", paste0(tmp, ".prune.in"),
                "--make-king-table", "--king-table-filter", "-1",
                "--out", tmp), "king")
    kin_f <- paste0(tmp, ".kin0")
    if (!file.exists(kin_f)) stop("plink2 produced no .kin0")

    k <- fread(kin_f); setnames(k, gsub("^#", "", names(k)))
    if (!nrow(k)) stop("plink2 .kin0 is empty: ", nrow(k), " pairs for ",
                       length(donors), " donors -- expected n*(n-1)/2")
    rel_max <- max(k$KINSHIP)
    message(sprintf("  KING: %d pairs, %d donors, max kinship %.4f",
                    nrow(k), length(donors), rel_max))

    k <- k[KINSHIP >= KIN_CUT]
    if (nrow(k)) {
        k[, donor1 := paste0(FID1, "::", IID1)][, donor2 := paste0(FID2, "::", IID2)]
        fl <- folds[, .(repeat_i, donor, outer_fold)]
        m1 <- merge(k, fl, by.x = "donor1", by.y = "donor", allow.cartesian = TRUE)
        m2 <- merge(m1, fl, by.x = c("donor2", "repeat_i"),
                    by.y = c("donor", "repeat_i"), suffixes = c("_1", "_2"))
        m2[, split := outer_fold_1 != outer_fold_2]
        rel_pairs <- m2[, .(donor1, donor2, KINSHIP, repeat_i,
                            fold1 = outer_fold_1, fold2 = outer_fold_2, split)]
        rel_summary <- rel_pairs[, .(n_repeats_split = sum(split), n_repeats = .N,
                                     kinship = KINSHIP[1]),
                                 by = .(donor1, donor2)][order(-kinship)]
    }
}

## --------------------------------------------------- C. cross-region concordance
## Same donors, same genome. A genuinely tagged locus is predictable everywhere.
concord <- data.table(); tail_ladder <- data.table()
ladder_one <- function(v, rid) data.table(
    run_id = rid, n = length(v),
    frac_gt_0.50 = mean(v > .50), frac_gt_0.70 = mean(v > .70),
    frac_gt_0.80 = mean(v > .80), frac_gt_0.85 = mean(v > .85),
    frac_gt_0.90 = mean(v > .90), frac_gt_0.95 = mean(v > .95),
    max = max(v))
tail_ladder <- ladder_one(dt$r2_pred_oof, RUN_ID)

if (length(COMPARE) > 0L) {
    hi_g <- GRanges(dt[r2_pred_oof > HI]$chrom,
                    IRanges(dt[r2_pred_oof > HI]$start + 1L, dt[r2_pred_oof > HI]$end))
    for (rid in COMPARE) {
        o <- read_vmrs(rid)
        tail_ladder <- rbind(tail_ladder, ladder_one(o$r2_pred_oof, rid))
        og <- GRanges(o$chrom, IRanges(o$start + 1L, o$end))

        ## the high tail specifically
        h <- as.data.table(findOverlaps(hi_g, og))
        h <- if (nrow(h)) h[, .(r2_other = max(o$r2_pred_oof[subjectHits])), by = queryHits] else
             data.table(queryHits = integer(), r2_other = numeric())
        ## and the whole catalog, for the rank correlation
        a <- as.data.table(findOverlaps(vmr, og))
        a <- a[, .(r2_other = max(o$r2_pred_oof[subjectHits])), by = queryHits]
        a[, r2_this := dt$r2_pred_oof[queryHits]]

        concord <- rbind(concord, data.table(
            other_run = rid,
            n_high_overlapping = nrow(h),
            n_high_total = length(hi_g),
            median_r2_other_for_high = if (nrow(h)) median(h$r2_other) else NA_real_,
            frac_other_gt_0.5_for_high = if (nrow(h)) mean(h$r2_other > .5) else NA_real_,
            frac_other_gt_0.9_for_high = if (nrow(h)) mean(h$r2_other > .9) else NA_real_,
            n_overlapping_all = nrow(a),
            spearman_all = cor(a$r2_this, a$r2_other, method = "spearman"),
            median_other_when_this_lt_0.05 = median(a[r2_this < .05]$r2_other)))
    }
}

## -------------------------------------------------------------------- output
out_dir <- file.path(ROOT, MODULE, "_m", "qc", RUN_ID)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(bins,  file.path(out_dir, "predictability-by-annotation-bin.tsv"), sep = "\t")
fwrite(tests, file.path(out_dir, "predictability-artifact-tests.tsv"), sep = "\t")
fwrite(mt,    file.path(out_dir, "predictability-length-matched.tsv"), sep = "\t")
fwrite(dt[, .(vmr_id, chrom, start, end, length_bp, r2_pred_oof, median_n_variants,
              segdup_frac, blacklist_frac, line_l1_frac, mappability)],
       file.path(out_dir, "per-locus-annotation.tsv.gz"), sep = "\t")
if (nrow(rel_summary)) {
    fwrite(rel_summary, file.path(out_dir, "relatedness-pairs.tsv"), sep = "\t")
    fwrite(rel_pairs, file.path(out_dir, "relatedness-fold-splits.tsv.gz"), sep = "\t")
}
fwrite(tail_ladder, file.path(out_dir, "r2-tail-ladder.tsv"), sep = "\t")
if (nrow(concord)) fwrite(concord, file.path(out_dir, "cross-region-concordance.tsv"), sep = "\t")

cat("\n=== r2_pred_oof bin x annotation (", RUN_ID, ") ===\n"); print(bins)
cat("\n=== logistic: high vs rest, raw and length/variant-adjusted ===\n"); print(tests)
cat("\n=== length-matched contrast (headline) ===\n"); print(mt)
cat("\n=== B. cryptic relatedness (KING >= ", KIN_CUT, ") ===\n", sep = "")
if (SKIP_REL) cat("  skipped (--skip-relatedness)\n") else
if (!nrow(rel_summary)) cat(sprintf(
    "  max observed kinship %.4f over all pairs, below the %.4f cutoff:\n  no relative pair exists for CV to split\n",
    rel_max, KIN_CUT)) else {
    print(rel_summary)
    cat("  pairs split across outer folds in >=1 repeat: ",
        rel_summary[n_repeats_split > 0, .N], " of ", nrow(rel_summary), "\n", sep = "")
}
cat("\n=== C. r2 tail ladder (threshold sensitivity) ===\n"); print(tail_ladder)
if (nrow(concord)) { cat("\n=== C. cross-region concordance ===\n"); print(concord) }
cat("\nWrote ", out_dir, "\n", sep = "")
