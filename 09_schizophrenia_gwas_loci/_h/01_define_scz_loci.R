#!/usr/bin/env Rscript
#### 09_schizophrenia_gwas_loci -- prespecify PGC3 loci (hg38) ####
##
## Usage:
##   Rscript _h/01_define_scz_loci.R --run-id scz-AA-caudate-YYYYMMDD
##
## AGENTS.md 7.8 makes locus definition a design invariant: PGC schizophrenia
## loci are defined independently of any methylation result. Nothing in this
## stage reads a VMR, a methylation value, or an upstream v2 run. It consumes
## only a published PGC3 locus table, the published index-SNP table, and the
## public wave-3 European summary statistics.
##
## Two locus sets (config locus_definition.locus_sets; the run manifest names
## the one used):
##   pgc3_published_gws   PRIMARY. PGC3 Supplementary Table 3: the 287 published
##                        genome-wide significant loci, as merged intervals,
##                        with their own index SNPs (282 autosomal).
##   pgc3_ld_clumps_p1e6  EXPLORATORY. PGC3 Supplementary Table 1: one interval
##                        per p < 1e-6 LD clump (612 autosomal), the set the
##                        scz-AA-*-20261008 runs used. Not loci, and not all
##                        genome-wide significant.
##
## The published index-SNP table and the summary statistics are hg19; the locus
## intervals, GTEx, and every v2 output are hg38. This stage is where that is
## reconciled, once, with liftOver, so no later stage has to guess a build.
## Variants that do not lift are DROPPED and COUNTED -- never silently carried
## with an hg19 coordinate into an hg38 join.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09_schizophrenia_gwas_loci", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
})

MODULE <- "09_schizophrenia_gwas_loci"
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

## ------------------------------------------------------------- locus set
## 00_new_run.R froze the locus set into the manifest; the config snapshot must
## agree, or the run would describe one set and analyse another.
manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
locus_set <- manifest$value[manifest$field == "locus_set"][1]
if (is.na(locus_set) || !nzchar(locus_set)) {
    stop("Run manifest records no locus_set; open the run with 00_new_run.R")
}
if (!identical(locus_set, ld$locus_set)) {
    stop("Manifest locus_set '", locus_set, "' differs from the config ",
         "snapshot's '", ld$locus_set, "'")
}
set_cfg <- ld$locus_sets[[locus_set]]
autosomes <- paste0("chr", 1:22)

## ------------------------------------------------------------------ index SNPs
## PGC3 Supplementary Table 1 carries the per-SNP statistics for both locus
## sets: 624 LD-clump index SNPs at p < 1e-6, hg19.
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
idx[, index_stats_source := "pgc3_st1_discovery"]

split_alleles <- function(a1a2) {
    a <- tstrsplit(as.character(a1a2), "/", fixed = TRUE)
    list(toupper(trimws(a[[1]])), toupper(trimws(a[[2]])))
}

if (identical(locus_set, "pgc3_published_gws")) {
    ## ---------------------------------------- published PGC3 loci (ST3)
    ## 287 merged regions, each containing >= 1 genome-wide significant SNP in
    ## the combined discovery-replication meta-analysis. The locus, its interval
    ## and its index SNPs are all PGC3's own: nothing is re-clumped or re-merged
    ## here, and each index SNP belongs to the locus that lists it.
    s3 <- as.data.table(readxl::read_excel(set_cfg$loci_xls, sheet = 1))
    need <- c("Chromosome", "top-index", "top-pos", "top-alleles", "top-P",
              "top-OR", "top-SE", "merge-LEFT", "merge-RIGHT", "indices", "genes")
    miss <- setdiff(need, names(s3))
    if (length(miss)) stop("Published locus table missing column(s): ",
                           paste(miss, collapse = ", "))
    if (nrow(s3) != as.integer(set_cfg$n_loci_published)) {
        stop("Published locus table has ", nrow(s3), " rows; config expects ",
             set_cfg$n_loci_published)
    }
    ## IDs follow the published row order (sorted by top-P), so they are
    ## fixed by PGC3 and cannot collide with the exploratory clump numbering.
    s3[, locus_id := sprintf("PGC3_%03d", .I)]
    s3[, chrom_hg19 := chr_norm(sub("^23$", "X", as.character(Chromosome)))]
    n_loci_source <- nrow(s3)
    n_loci_nonautosomal <- sum(!s3$chrom_hg19 %in% autosomes)
    ## Autosomes only: PGC3 wave-3 public sumstats are autosomal, and Module 05
    ## mapped autosomes only (config/meqtl_parameters.yml:chromosomes: autosomal).
    s3 <- s3[chrom_hg19 %in% autosomes]

    ## Lift both interval ends. An end that falls in sequence deleted from hg38
    ## (3 of the 564 autosomal ends: two beside the chr5 centromere, one at the
    ## chr16 telomere) is walked INWARD in 1 kb steps to the first position that
    ## lifts, so the interval shrinks at its unmappable edge rather than the
    ## published locus being lost. The shift is recorded per locus. A locus
    ## whose end cannot be placed within max_end_walk_bp is dropped and counted.
    step_bp <- 1000L
    max_end_walk_bp <- 100000L
    ends <- rbind(
        s3[, .(locus_id, side = "L", chrom_hg19, pos0 = as.integer(`merge-LEFT`),
               inward = 1L)],
        s3[, .(locus_id, side = "R", chrom_hg19, pos0 = as.integer(`merge-RIGHT`),
               inward = -1L)])
    cand <- ends[, .(offset = seq(0L, max_end_walk_bp, by = step_bp)),
                 by = .(locus_id, side, chrom_hg19, pos0, inward)]
    cand[, `:=`(pos_hg19 = pos0 + inward * offset,
                end_key = paste0(locus_id, ":", side, ":", offset))]
    ends_lift <- lift_to_hg38(cand, "end_key", file.path(run_dir, "gwas"), "locus_ends")
    setorder(ends_lift, locus_id, side, offset)
    ends_lift <- unique(ends_lift, by = c("locus_id", "side"))
    iv <- dcast(ends_lift, locus_id + chrom_hg38 ~ side,
                value.var = c("pos_hg38", "offset"))
    iv <- iv[is.finite(pos_hg38_L) & is.finite(pos_hg38_R)]
    iv[, `:=`(start = pmin(pos_hg38_L, pos_hg38_R), end = pmax(pos_hg38_L, pos_hg38_R),
              end_walked_inward_bp = offset_L + offset_R)]
    loci <- merge(s3, iv[, .(locus_id, chrom = chrom_hg38, start, end,
                             end_walked_inward_bp)],
                  by = "locus_id")
    n_loci_lift_dropped <- nrow(s3) - nrow(loci)
    n_loci_end_walked <- sum(loci$end_walked_inward_bp > 0)
    loci[, `:=`(span_hg19 = `merge-RIGHT` - `merge-LEFT`, span_hg38 = end - start)]

    ## Index SNPs, one row per (locus, SNP), statistics from ST1. The one
    ## published index SNP absent from ST1 is a top index SNP, so ST3's own
    ## top-* columns describe it.
    pub <- loci[, .(index_snp = trimws(unlist(strsplit(as.character(indices), ",")))),
                by = locus_id]
    pub <- pub[nzchar(index_snp)]
    if (anyDuplicated(pub$index_snp)) {
        stop("An index SNP is listed under more than one published locus: ",
             paste(head(pub$index_snp[duplicated(pub$index_snp)]), collapse = ", "))
    }
    st1 <- idx[, .(index_snp, chrom_hg19 = chr_norm(chrom_hg19), pos_hg19,
                   gwas_p, gwas_or, gwas_se,
                   risk_allele, other_allele, index_stats_source)]
    pub <- merge(pub, st1, by = "index_snp", all.x = TRUE)
    top <- loci[, .(locus_id, index_snp = `top-index`, chrom_hg19,
                    pos_hg19 = as.integer(`top-pos`), gwas_p = as.numeric(`top-P`),
                    gwas_or = as.numeric(`top-OR`), gwas_se = as.numeric(`top-SE`),
                    a1a2 = `top-alleles`)]
    top[, c("risk_allele", "other_allele") := split_alleles(a1a2)]
    no_st1 <- pub[is.na(pos_hg19), index_snp]
    if (length(setdiff(no_st1, top$index_snp))) {
        stop("Published index SNP(s) with no statistics in ST1 or ST3: ",
             paste(setdiff(no_st1, top$index_snp), collapse = ", "))
    }
    fill <- top[index_snp %in% no_st1, !"a1a2"][, index_stats_source := "pgc3_st3_combined"]
    pub <- rbind(pub[!is.na(pos_hg19)], fill, use.names = TRUE)

    idx_lift <- lift_to_hg38(pub, "index_snp", file.path(run_dir, "gwas"), "index_snps")
    n_index_in <- attr(idx_lift, "n_input")
    index_out <- idx_lift[, .(locus_id, index_snp,
                              chrom = chrom_hg38, pos_hg38,
                              risk_allele, other_allele,
                              gwas_p, gwas_or, gwas_se,
                              chrom_hg19 = chr_norm(chrom_hg19), pos_hg19,
                              index_stats_source)]
    lk <- loci[, .(locus_id, l_chrom = chrom, l_start = start, l_end = end)]
    index_out <- merge(index_out, lk, by = "locus_id", all.x = TRUE)
    index_out[, index_within_interval := chrom == l_chrom & pos_hg38 >= l_start &
                                         pos_hg38 <= l_end]
    index_out[, c("l_chrom", "l_start", "l_end") := NULL]

    ## The locus's top index SNP and its PUBLISHED combined p are the locus-level
    ## values (prioritization rank key 3).
    top_pos <- index_out[, .(index_snp, index_pos_hg38 = pos_hg38)]
    loci_out <- loci[, .(locus_id, chrom, start, end,
                         index_snp = `top-index`,
                         pgc3_index_pvalue = as.numeric(`top-P`),
                         pgc3_index_or = as.numeric(`top-OR`),
                         n_index_snps = lengths(strsplit(as.character(indices), ",")),
                         pgc3_genes = as.character(genes),
                         span_hg19, span_hg38, end_walked_inward_bp)]
    loci_out <- merge(loci_out, top_pos, by = "index_snp", all.x = TRUE)
    loci_sources <- set_cfg$loci_xls
} else if (identical(locus_set, "pgc3_ld_clumps_p1e6")) {
    ## ------------------------------- exploratory: ST1 LD clumps, p < 1e-6
    ## One interval per clump, as the accepted scz-AA-*-20261008 runs used them.
    loci <- fread(set_cfg$pgc3_loci_hg38)
    for (c_ in c("chrom", "start", "end", "locus")) {
        if (!c_ %in% names(loci)) stop("loci table missing column ", c_)
    }
    loci[, chrom := chr_norm(chrom)]
    loci[, locus_id := as.character(locus)]
    loci[, `:=`(start = as.integer(start), end = as.integer(end))]
    n_loci_source <- nrow(loci)
    n_loci_nonautosomal <- sum(!loci$chrom %in% autosomes)
    n_loci_lift_dropped <- 0L
    n_loci_end_walked <- 0L
    loci <- loci[chrom %in% autosomes]
    setorder(loci, chrom, start)

    idx_lift <- lift_to_hg38(idx, "index_snp", file.path(run_dir, "gwas"), "index_snps")
    n_index_in <- attr(idx_lift, "n_input")

    ## Assign each index SNP to the smallest interval containing it. An index
    ## SNP outside every interval keeps locus_id NA rather than being forced
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
                              chrom_hg19 = chr_norm(chrom_hg19), pos_hg19,
                              index_stats_source)]
    index_out[, index_within_interval := !is.na(locus_id)]

    ## One index SNP per locus for the prioritization tie-break: the most
    ## significant, chosen by p-value alone and therefore independent of
    ## anything this module computes.
    best_idx <- index_out[!is.na(locus_id) & is.finite(gwas_p)]
    setorder(best_idx, locus_id, gwas_p)
    best_idx <- unique(best_idx, by = "locus_id")
    loci_out <- merge(loci[, .(locus_id, chrom, start, end)],
                      best_idx[, .(locus_id, index_snp,
                                   index_pos_hg38 = pos_hg38,
                                   pgc3_index_pvalue = gwas_p,
                                   pgc3_index_or = gwas_or)],
                      by = "locus_id", all.x = TRUE)
    loci_sources <- set_cfg$pgc3_loci_hg38
} else {
    stop("No stage 01 reader for locus set '", locus_set, "'")
}

setorder(index_out, chrom, pos_hg38)
write_atomic(index_out, file.path(run_dir, "results", "scz-index-snps.tsv"))

loci_out[, `:=`(window_start = pmax(1L, start - flank), window_end = end + flank,
                locus_set = locus_set)]
setcolorder(loci_out, c("locus_id", "chrom", "start", "end", "window_start",
                        "window_end", "index_snp", "index_pos_hg38",
                        "pgc3_index_pvalue", "pgc3_index_or"))
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
## n_loci_published is the number of loci that ENTER the analysis (autosomal,
## lifted); the gate's min_loci_tested reads it. The source-table count and
## what was dropped on the way are recorded beside it.
summary_dt <- data.table(
    run_id = opts$run_id,
    locus_set = locus_set,
    locus_set_role = set_cfg$role,
    n_loci_in_source_table = n_loci_source,
    n_loci_nonautosomal_dropped = n_loci_nonautosomal,
    n_loci_lift_dropped = n_loci_lift_dropped,
    n_loci_end_walked_inward = n_loci_end_walked,
    n_loci_published = nrow(loci_out),
    n_loci_with_index_snp = sum(!is.na(loci_out$index_snp)),
    n_index_snps_published = n_index_in,
    n_index_snps_lifted = nrow(idx_lift),
    n_index_snps_assigned_to_locus = sum(!is.na(index_out$locus_id)),
    n_index_snps_within_interval = sum(index_out$index_within_interval, na.rm = TRUE),
    n_gwas_variants_input = n_gw_in,
    n_gwas_variants_lifted_into_windows = n_gw_lifted,
    n_gwas_variants_in_locus_windows = nrow(gw_win),
    locus_flank_bp = flank,
    genome_build = "hg38",
    loci_source = loci_sources,
    index_source = ld$pgc3_index_snp_xls,
    sumstats_source = ld$pgc3_sumstats_hg19,
    methylation_used = FALSE
)
write_atomic(summary_dt, file.path(run_dir, "results", "locus-definition-summary.tsv"))
print(summary_dt[, .(n_loci_published, n_loci_with_index_snp,
                     n_gwas_variants_in_locus_windows)])
message("[09] loci defined (", locus_set, "): ", nrow(loci_out), " intervals, ",
        nrow(gw_win), " GWAS variants in window")
