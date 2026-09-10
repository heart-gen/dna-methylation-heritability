#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- prespecify PGC3 loci (hg38) ####
##
## Usage:
##   Rscript _h/01_define_scz_loci.R --run-id scz-AA-caudate-YYYYMMDD
##
## AGENTS.md 7.8 makes locus definition a design invariant: PGC schizophrenia
## loci are defined independently of any methylation result. Nothing in this
## stage reads a VMR, a methylation value, or an upstream v2 run. It consumes
## only the published fine-mapped locus intervals, the published index-SNP
## table, and the public wave-3 European summary statistics.
##
## The published index-SNP table and the summary statistics are hg19; the locus
## intervals, GTEx, and every v2 output are hg38. This stage is where that is
## reconciled, once, with liftOver, so no later stage has to guess a build.
## Variants that do not lift are DROPPED and COUNTED -- never silently carried
## with an hg19 coordinate into an hg38 join.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09_schizophrenia_risk_application", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
})

MODULE <- "09_schizophrenia_risk_application"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

scz <- load_run_config("schizophrenia", run_dir)
ld <- scz$locus_definition
flank <- as.integer(ld$locus_flank_bp)

liftover_bin <- Sys.getenv("V2_LIFTOVER",
                           "/projects/p32505/opt/envs/genomics/bin/liftOver")
if (!file.exists(liftover_bin)) stop("liftOver not found: ", liftover_bin)

`%||%` <- function(a, b) if (is.null(a)) b else a
chr_norm <- function(x) paste0("chr", sub("^chr", "", as.character(x)))

## ------------------------------------------------------------ liftOver hg19
## Lift a table carrying chrom_hg19/pos_hg19 and a unique `key`. Returns the
## table with pos_hg38 added, keeping only rows that lifted to the SAME
## chromosome -- a cross-chromosome lift is a mapping artefact, not a coordinate.
lift_to_hg38 <- function(dt, key_col, work_dir, tag) {
    dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
    bed_in <- file.path(work_dir, paste0(tag, ".hg19.bed"))
    bed_out <- file.path(work_dir, paste0(tag, ".hg38.bed"))
    unmapped <- file.path(work_dir, paste0(tag, ".unmapped.bed"))

    src <- dt[is.finite(pos_hg19) & !is.na(chrom_hg19)]
    fwrite(src[, .(chr_norm(chrom_hg19), pos_hg19 - 1L, pos_hg19,
                   get(key_col))],
           bed_in, sep = "\t", col.names = FALSE, quote = FALSE)
    st <- system2(liftover_bin, c(shQuote(bed_in),
                                  shQuote(ld$liftover_chain_hg19_to_hg38),
                                  shQuote(bed_out), shQuote(unmapped)),
                  stdout = FALSE, stderr = FALSE)
    if (st != 0 || !file.exists(bed_out)) {
        stop("liftOver failed for ", tag, " (exit ", st, ")")
    }
    lifted <- fread(bed_out, header = FALSE,
                    col.names = c("chrom_hg38", "start0", "pos_hg38", key_col))
    if (nrow(lifted) == 0) stop("liftOver mapped no ", tag, " records")
    lifted <- unique(lifted, by = key_col)
    out <- merge(dt, lifted[, c(key_col, "chrom_hg38", "pos_hg38"),
                            with = FALSE],
                 by = key_col, all.x = FALSE)
    same_chr <- chr_norm(out$chrom_hg19) == out$chrom_hg38
    dropped_chr <- sum(!same_chr)
    out <- out[same_chr]
    attr(out, "n_input") <- nrow(src)
    attr(out, "n_lifted") <- nrow(out)
    attr(out, "n_cross_chrom_dropped") <- dropped_chr
    out[]
}

## ------------------------------------------------------------- locus intervals
loci <- fread(ld$pgc3_loci_hg38)
for (c_ in c("chrom", "start", "end", "locus")) {
    if (!c_ %in% names(loci)) stop("loci table missing column ", c_)
}
loci[, chrom := chr_norm(chrom)]
loci[, locus_id := as.character(locus)]
loci[, `:=`(start = as.integer(start), end = as.integer(end))]
## Autosomes only: PGC3 wave-3 public sumstats are autosomal, and Module 05
## mapped autosomes only (config/meqtl_parameters.yml:chromosomes: autosomal).
loci <- loci[chrom %in% paste0("chr", 1:22)]
loci[, `:=`(window_start = pmax(1L, start - flank), window_end = end + flank)]
setorder(loci, chrom, start)

## ------------------------------------------------------------------ index SNPs
idx <- as.data.table(readxl::read_excel(ld$pgc3_index_snp_xls))
index_cols <- c(SNP = "index_snp", CHR = "chrom_hg19", BP = "pos_hg19",
                P = "gwas_p", OR = "gwas_or", SE = "gwas_se", A1A2 = "a1a2")
missing_idx <- setdiff(names(index_cols), names(idx))
if (length(missing_idx)) {
    stop("Index-SNP table missing column(s): ", paste(missing_idx, collapse = ", "))
}
setnames(idx, names(index_cols), unname(index_cols))
idx[, `:=`(index_snp = as.character(index_snp),
           pos_hg19 = as.integer(pos_hg19),
           gwas_p = as.numeric(gwas_p),
           gwas_or = as.numeric(gwas_or),
           gwas_se = as.numeric(gwas_se))]
alleles <- tstrsplit(as.character(idx$a1a2), "/", fixed = TRUE)
idx[, `:=`(risk_allele = toupper(trimws(alleles[[1]])),
           other_allele = toupper(trimws(alleles[[2]])))]
idx <- unique(idx, by = "index_snp")

idx_lift <- lift_to_hg38(idx, "index_snp", file.path(run_dir, "gwas"), "index_snps")
n_index_in <- attr(idx_lift, "n_input")

## Assign each index SNP to the smallest published interval containing it. An
## index SNP outside every interval keeps locus_id NA rather than being forced
## into the nearest locus.
setkey(loci, chrom, start, end)
idx_lift[, `:=`(chrom = chrom_hg38, s = pos_hg38, e = pos_hg38)]
hit <- foverlaps(idx_lift[, .(index_snp, chrom, start = s, end = e)],
                 loci[, .(chrom, start, end, locus_id, span = end - start)],
                 by.x = c("chrom", "start", "end"), type = "within",
                 nomatch = NULL)
setorder(hit, index_snp, span)
hit <- unique(hit, by = "index_snp")
idx_lift <- merge(idx_lift, hit[, .(index_snp, locus_id)],
                  by = "index_snp", all.x = TRUE)

index_out <- idx_lift[, .(locus_id, index_snp,
                          chrom = chrom_hg38, pos_hg38,
                          risk_allele, other_allele,
                          gwas_p, gwas_or, gwas_se,
                          chrom_hg19 = chr_norm(chrom_hg19), pos_hg19)]
setorder(index_out, chrom, pos_hg38)
write_atomic(index_out, file.path(run_dir, "results", "scz-index-snps.tsv"))

## One index SNP per locus for the prioritization tie-break: the most
## significant, chosen by p-value alone and therefore independent of anything
## this module computes.
best_idx <- index_out[!is.na(locus_id) & is.finite(gwas_p)]
setorder(best_idx, locus_id, gwas_p)
best_idx <- unique(best_idx, by = "locus_id")
loci_out <- merge(loci[, .(locus_id, chrom, start, end, window_start, window_end)],
                  best_idx[, .(locus_id, index_snp,
                               index_pos_hg38 = pos_hg38,
                               pgc3_index_pvalue = gwas_p,
                               pgc3_index_or = gwas_or)],
                  by = "locus_id", all.x = TRUE)
setorder(loci_out, chrom, start)
write_atomic(loci_out, file.path(run_dir, "results", "scz-loci.tsv"))

## ------------------------------------------------------- GWAS summary stats
## Sliced to the locus windows before anything else touches it: coloc only ever
## needs the regions, and the whole-genome file is 7.6M variants.
##
## This is done as a shell pipeline rather than in R on purpose. Reading the
## full table, lifting it, and joining in memory was OOM-killed on the submit
## host, and this stage runs on the submit host by design (every array task
## reads its output). awk streams the coordinates out, liftOver maps them,
## bedtools selects the windows, and only the surviving rows are ever read into
## R. Line numbers of the de-commented stream are the join key, so no
## coordinate or allele string has to round-trip through two formats.
gwas_dir <- file.path(run_dir, "gwas")
win_bed <- file.path(gwas_dir, "locus-windows.bed")
fwrite(loci_out[, .(chrom, window_start - 1L, window_end, locus_id)],
       win_bed, sep = "\t", col.names = FALSE, quote = FALSE)

bedtools_bin <- Sys.getenv("V2_BEDTOOLS",
                           "/projects/p32505/opt/envs/genomics/bin/bedtools")
if (!file.exists(bedtools_bin)) stop("bedtools not found: ", bedtools_bin)

slicer <- file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(repo_root(), MODULE, "_h")), "slice_gwas_sumstats.sh")
if (!file.exists(slicer)) stop("Missing slicing script: ", slicer)
st <- system2("bash", shQuote(c(slicer, ld$pgc3_sumstats_hg19,
                                ld$liftover_chain_hg19_to_hg38, win_bed,
                                gwas_dir, liftover_bin, bedtools_bin)))
if (st != 0) stop("GWAS slicing pipeline failed (exit ", st, ")")

## hg38 coordinate and locus assignment for each surviving row, keyed by the
## same line number.
lifted <- fread(file.path(gwas_dir, "sumstats.inwindow.bed"), header = FALSE,
                col.names = c("chrom_hg38", "start0", "pos_hg38", "row_key",
                              "win_chrom", "win_start", "win_end", "locus_id"),
                colClasses = list(character = c(1, 5, 8)))
n_gw_lifted <- nrow(unique(lifted, by = "row_key"))

gw <- fread(file.path(gwas_dir, "sumstats.selected.tsv"), showProgress = FALSE)
req <- c("CHROM", "ID", "POS", "A1", "A2", "BETA", "SE", "PVAL", "NEFF")
missing <- setdiff(req, names(gw))
if (length(missing)) stop("GWAS sumstats missing column(s): ",
                          paste(missing, collapse = ", "))
if (!"row_key" %in% names(gw)) {
    stop("sumstats.selected.tsv has no row_key column; the slicing script must ",
         "emit it (the selected file is a subset, so row order is not a key).")
}
gw <- gw[, .(row_key = as.integer(row_key), variant_rsid = as.character(ID),
             chrom_hg19 = chr_norm(CHROM), pos_hg19 = as.integer(POS),
             effect_allele = toupper(A1), other_allele = toupper(A2),
             beta = as.numeric(BETA), se = as.numeric(SE),
             pvalue = as.numeric(PVAL), neff = as.numeric(NEFF))]
gw <- gw[is.finite(beta) & is.finite(se) & se > 0]

gw_win <- merge(gw, lifted[, .(row_key, chrom = chrom_hg38, pos_hg38, locus_id)],
                by = "row_key", allow.cartesian = TRUE)
## A lift that changed chromosome is a mapping artefact, not a coordinate.
gw_win <- gw_win[chrom == chrom_hg19]
gw_win[, variant_key := paste0(chrom, "_", pos_hg38)]
gw_win <- unique(gw_win, by = c("locus_id", "variant_key",
                                "effect_allele", "other_allele"))
setorder(gw_win, chrom, pos_hg38, locus_id)
fwrite(gw_win[, .(variant_rsid, chrom, pos_hg38, effect_allele, other_allele,
                  beta, se, pvalue, neff, locus_id, variant_key)],
       file.path(gwas_dir, "gwas-sumstats-hg38.tsv.gz"),
       sep = "\t", compress = "gzip")
n_gw_in <- as.integer(strsplit(trimws(system2(
    "wc", c("-l", shQuote(file.path(gwas_dir, "sumstats.hg19.bed"))),
    stdout = TRUE)), " +")[[1]][1])

## ------------------------------------------------------------------ summary
summary_dt <- data.table(
    run_id = opts$run_id,
    n_loci_published = nrow(loci),
    n_loci_with_index_snp = sum(!is.na(loci_out$index_snp)),
    n_index_snps_published = n_index_in,
    n_index_snps_lifted = nrow(idx_lift),
    n_index_snps_assigned_to_locus = sum(!is.na(index_out$locus_id)),
    n_gwas_variants_input = n_gw_in,
    n_gwas_variants_lifted_into_windows = n_gw_lifted,
    n_gwas_variants_in_locus_windows = nrow(gw_win),
    locus_flank_bp = flank,
    genome_build = "hg38",
    loci_source = ld$pgc3_loci_hg38,
    index_source = ld$pgc3_index_snp_xls,
    sumstats_source = ld$pgc3_sumstats_hg19,
    methylation_used = FALSE
)
write_atomic(summary_dt, file.path(run_dir, "results", "locus-definition-summary.tsv"))
print(summary_dt[, .(n_loci_published, n_loci_with_index_snp,
                     n_gwas_variants_in_locus_windows)])
message("[09] loci defined: ", nrow(loci), " published intervals, ",
        nrow(gw_win), " GWAS variants in window")
