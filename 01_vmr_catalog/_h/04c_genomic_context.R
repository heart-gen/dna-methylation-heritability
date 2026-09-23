#### 01_vmr_catalog / 04c_genomic_context: where VMRs sit in the genome ####
##
## Descriptive characterization of the corrected catalog: which genic
## compartment each VMR overlaps, and how far it is from the nearest gene.
##
## This is NOT an enrichment analysis. No comparison group, no test, no model.
## The primary genomic-enrichment model is a locked PI decision (AGENTS.md 12)
## and repeat/repressive architecture is Module 04's remit; this stage only
## describes the catalog, the way width and CpG count do.
##
## The legacy equivalent
## (local-snp-prediction/BA_only/tissue_comparison/annotation/_h/01.annotate_vmrs.R)
## cannot be reused: it reads vmr-analysis/{tissue}/_m/vmr/chr_*/vmr.bed, the
## catalog invalidated by the V1 donor-row misalignment (AGENTS.md 8).
##
## Compartments are assigned by priority, not by multiple membership, so the
## proportions sum to 1 and the panel is readable: a VMR overlapping both a
## promoter and an intron is called a promoter.
##
## F14: the priority order is now separate from the display order, because they
## are not the same thing and conflating them made one category impossible.
## `Exon` used to be tested before `3' UTR`. A 3' UTR *is* an exonic interval, so
## every VMR that could have been called `3' UTR` was called `Exon` first and the
## category could only ever be zero -- which then reached Figure 1 as an empty
## bar, asserting "we looked and found none" where the truth was "we could not
## look". Both UTRs now precede the generic exon layer, which is the order the
## 5' UTR already had, and `assert_levels_reachable()` below fails the stage if
## any declared category is unreachable again. `Exon` consequently means an
## exonic interval that is not a UTR, mostly coding sequence.
##
## Usage:
##   Rscript 04c_genomic_context.R --cohort AA --region caudate --run-id ID

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

opts <- parse_v2_args(require = c("cohort", "region", "run_id"))
cohort <- opts$cohort; region <- opts$region

module_root <- file.path(V2_ROOT, "01_vmr_catalog")
run_dir <- file.path(module_root, "_m", "runs", opts$run_id)
qc_dir  <- file.path(run_dir, "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

vmr <- fread(file.path(run_dir, "vmr", "vmr_catalog.tsv"))
vmr_set_id <- vmr$vmr_set_id[1]
g_vmr <- GRanges(vmr$chr, IRanges(vmr$start, vmr$end))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
keep <- paste0("chr", chrom_order(include_sex = FALSE))
trim_to <- function(g) {
    g <- g[as.character(seqnames(g)) %in% keep]
    GenomicRanges::reduce(GenomicRanges::granges(g), ignore.strand = TRUE)
}

genes_gr <- suppressMessages(GenomicFeatures::genes(txdb))
genes_gr <- genes_gr[as.character(seqnames(genes_gr)) %in% keep]

## Promoter defined as 2 kb upstream to 200 bp downstream of the TSS, the
## annotatr hg38_basicgenes convention the legacy script used.
prom  <- trim_to(GenomicFeatures::promoters(txdb, upstream = 2000, downstream = 200))
utr5  <- trim_to(unlist(GenomicFeatures::fiveUTRsByTranscript(txdb)))
utr3  <- trim_to(unlist(GenomicFeatures::threeUTRsByTranscript(txdb)))
exons <- trim_to(GenomicFeatures::exons(txdb))
introns <- trim_to(unlist(GenomicFeatures::intronsByTranscript(txdb)))

## Display order, 5' to 3' along a gene. This is what the figure axis uses; it
## says nothing about which category wins an overlap.
LEVELS <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Intergenic")

## Assignment order, most specific first. `Intergenic` is the fallthrough and is
## therefore not a layer.
PRIORITY <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron")
layers <- list(Promoter = prom, `5' UTR` = utr5, `3' UTR` = utr3,
               Exon = exons, Intron = introns)
stopifnot(identical(sort(c(PRIORITY, "Intergenic")), sort(LEVELS)),
          identical(names(layers), PRIORITY))

#' Fail if a declared category cannot be assigned to any interval in the genome.
#'
#' A priority classification can silently swallow a category whenever an earlier
#' layer contains a later one. The check is on the annotation, not on this
#' catalog's VMRs: a category that no VMR happens to hit is a real zero, while a
#' category no interval anywhere could hit is a coding defect, and a count of
#' zero looks identical either way.
assert_levels_reachable <- function(layers, fallthrough, seq_keep) {
    covered <- NULL
    for (nm in names(layers)) {
        unique_bp <- if (is.null(covered)) {
            sum(width(layers[[nm]]))
        } else {
            sum(width(GenomicRanges::setdiff(layers[[nm]], covered,
                                             ignore.strand = TRUE)))
        }
        if (unique_bp == 0) {
            stop("Genomic-context category '", nm, "' is unreachable: every ",
                 "interval it could claim is already claimed by an earlier ",
                 "layer (", paste(names(layers)[seq_len(match(nm, names(layers)) - 1L)],
                                  collapse = ", "), ").\n  A category that can ",
                 "only ever be zero must be reordered or removed, never plotted.")
        }
        message("[reach] ", nm, ": ", format(unique_bp, big.mark = ","),
                " bp not claimed by a higher-priority layer")
        covered <- if (is.null(covered)) layers[[nm]] else
            GenomicRanges::reduce(c(covered, layers[[nm]]), ignore.strand = TRUE)
    }
    ## The fallthrough is reachable iff the layers do not tile the chromosomes.
    genome_bp <- sum(as.numeric(seqlengths(txdb)[seq_keep]), na.rm = TRUE)
    if (genome_bp - sum(as.numeric(width(covered))) <= 0) {
        stop("Genomic-context category '", fallthrough, "' is unreachable: the ",
             "annotation layers cover every base of the analyzed chromosomes.")
    }
    invisible(TRUE)
}
assert_levels_reachable(layers, "Intergenic", keep)

assigned <- rep(NA_character_, length(g_vmr))
for (nm in PRIORITY) {
    hit <- is.na(assigned) & IRanges::overlapsAny(g_vmr, layers[[nm]])
    assigned[hit] <- nm
}
assigned[is.na(assigned)] <- "Intergenic"

## Every declared level gets a row, including a true zero. The old table emitted
## only observed levels, so an absent category and a category with no VMRs were
## indistinguishable downstream -- which is how the unreachable one went unnoticed.
observed <- data.table(genomic_context = factor(assigned, levels = LEVELS))[
    , .(n = .N), by = genomic_context]
context <- data.table(
    cohort = cohort, region = region, vmr_set_id = vmr_set_id,
    genomic_context = factor(LEVELS, levels = LEVELS), n_vmrs = 0L)
context[observed, n_vmrs := i.n, on = "genomic_context"]
context[, frac_vmrs := n_vmrs / sum(n_vmrs)]
setorder(context, genomic_context)
stopifnot(abs(sum(context$frac_vmrs) - 1) < 1e-9,
          nrow(context) == length(LEVELS))

## Distance to the nearest gene, per VMR, for the distribution panel.
d <- distanceToNearest(g_vmr, genes_gr, ignore.strand = TRUE)
dist_dt <- data.table(
    cohort = cohort, region = region, vmr_set_id = vmr_set_id,
    vmr_id = vmr$vmr_id[queryHits(d)],
    genomic_context = assigned[queryHits(d)],
    distance_to_nearest_gene = mcols(d)$distance)

write_atomic(context, file.path(qc_dir, "genomic_context.tsv"))
write_atomic(dist_dt, file.path(qc_dir, "distance_to_nearest_gene.tsv"))

print(context[, .(genomic_context, n_vmrs, frac = round(frac_vmrs, 4))])
message("[done] genomic context written to ", qc_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
