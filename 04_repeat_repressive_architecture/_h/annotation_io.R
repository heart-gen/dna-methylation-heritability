## Annotation readers for 04_repeat_repressive_architecture.
##
## Two builds are in play and conflating them would silently destroy the result.
## VMRs are hg38. The Roadmap consolidated epigenomes -- the only brain-specific
## H3K9me3 and ChromHMM tracks the project has -- are hg19. Every function here
## states which build it returns, and `lift_vmrs_to_hg19()` is the single place
## the crossing happens.

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(rtracklayer)
})

#' Resolve a config path that may be repo-relative or absolute.
annot_path <- function(p, root = repo_root()) {
    if (is.null(p) || !nzchar(p)) stop("Missing annotation path in config")
    out <- if (startsWith(p, "/")) p else file.path(root, p)
    if (!file.exists(out)) stop("Annotation file not found: ", out)
    out
}

#' Roadmap epigenome ID for a brain region (hg19 tracks).
roadmap_eid <- function(annot, region) {
    eid <- annot$chromatin$roadmap_eid[[region]]
    if (is.null(eid)) {
        stop("config/repeat_annotations.yml has no Roadmap EID for region '",
             region, "'. It is a PI-lock key.")
    }
    eid
}

#' A gappedPeak chromatin track for one epigenome. Returns hg19 GRanges.
#'
#' gappedPeak is BED12-like: 15 columns, with the peak span in 1-3 and a block
#' structure in 10-12. rtracklayer::import() will not parse it as BED, so this
#' reads columns 1-3 with fread, exactly as the legacy analysis did. The block
#' structure is deliberately unused, so the "overlap" being tested is overlap
#' with the called peak, not with its sub-blocks.
#'
#' `key` names the block under `chromatin:` in config (h3k9me3, h3k27me3,
#' h3k27ac). Parameterized rather than one function per mark so that adding a
#' track is a config edit, not a code edit.
load_gappedpeak_hg19 <- function(annot, region, key) {
    eid <- roadmap_eid(annot, region)
    cfg <- annot$chromatin[[key]]
    if (is.null(cfg)) stop("No chromatin track configured under key: ", key)
    f <- annot_path(file.path(cfg$dir, gsub("{eid}", eid, cfg$template, fixed = TRUE)))
    dt <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE, select = 1:3,
                col.names = c("chrom", "start", "end"))
    if (nrow(dt) == 0) stop(key, " track is empty: ", f)
    GRanges(dt$chrom, IRanges(dt$start + 1L, dt$end))   # BED is 0-based half-open
}

#' A set of 15-state ChromHMM mnemonic states for one epigenome. Returns hg19
#' GRanges.
#'
#' `key` names the block under `chromatin:` in config (quiescent, bivalent,
#' accessible); the states themselves are listed there, never hardcoded here.
#' All three read the same mnemonics file and differ only in the state filter.
load_chromhmm_states_hg19 <- function(annot, region, key) {
    eid <- roadmap_eid(annot, region)
    cfg <- annot$chromatin[[key]]
    if (is.null(cfg)) stop("No chromatin track configured under key: ", key)
    f <- annot_path(file.path(cfg$dir, gsub("{eid}", eid, cfg$template, fixed = TRUE)))
    dt <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE,
                col.names = c("chrom", "start", "end", "state"))
    states <- unlist(cfg$states)
    keep <- dt[state %in% states]
    ## A typo in a state name would otherwise produce an empty annotation and a
    ## silent null result, which is indistinguishable from a real null.
    missing <- setdiff(states, unique(dt$state))
    if (length(missing) > 0) {
        stop("ChromHMM state(s) not present in ", f, ": ",
             paste(missing, collapse = ", "), ". Observed states: ",
             paste(sort(unique(dt$state)), collapse = ", "))
    }
    if (nrow(keep) == 0) {
        stop("No ChromHMM segment matched state(s) ", paste(states, collapse = ","),
             " in ", f)
    }
    GRanges(keep$chrom, IRanges(keep$start + 1L, keep$end))
}

#' Lift hg38 VMRs down to hg19.
#'
#' Returns a list with `gr` (hg19 GRanges carrying vmr_id) and `report`, a row
#' per VMR recording how many hg19 intervals it produced. A VMR that fails to
#' lift, or lifts to more than one interval, is DROPPED and reported -- never
#' silently collapsed to its first interval, which would quietly relocate a
#' locus and then test it against the wrong chromatin.
lift_vmrs_to_hg19 <- function(vmr_hg38, annot) {
    chain_f <- annot_path(annot$chromatin$liftover$chain)
    chain <- import.chain(chain_f)
    lifted <- liftOver(vmr_hg38, chain)
    n_int <- elementNROWS(lifted)

    report <- data.table(
        vmr_id = mcols(vmr_hg38)$vmr_id,
        n_hg19_intervals = as.integer(n_int),
        liftover_status = fifelse(n_int == 0L, "unmapped",
                           fifelse(n_int == 1L, "unique", "multi_mapping"))
    )
    keep <- n_int == 1L
    gr <- unlist(lifted[keep])
    mcols(gr)$vmr_id <- mcols(vmr_hg38)$vmr_id[keep]
    ## Width can change slightly across builds; overlap fractions must be taken
    ## against the hg19 width, which is what actually intersected the track.
    list(gr = gr, report = report)
}

#' Summed overlap fraction and any-overlap flag, per query interval.
#'
#' Both are kept: the fraction is the more informative outcome but is sensitive
#' to interval length (hence the length covariate), while the binary call is
#' what the legacy analysis reported and is needed for old-vs-new comparison.
overlap_features <- function(gr, ann_gr, prefix) {
    ann <- GenomicRanges::reduce(ann_gr, ignore.strand = TRUE)
    hits <- findOverlaps(gr, ann, ignore.strand = TRUE)
    frac <- numeric(length(gr))
    if (length(hits) > 0) {
        inter <- pintersect(gr[queryHits(hits)], ann[subjectHits(hits)],
                            ignore.strand = TRUE)
        cov <- tapply(width(inter), queryHits(hits), sum)
        frac[as.integer(names(cov))] <- as.numeric(cov)
    }
    frac <- pmin(frac / width(gr), 1)
    out <- data.table(frac, frac > 0)
    setnames(out, c(paste0(prefix, "_frac"), paste0(prefix, "_any")))
    out
}

#' Overlap features from a BED-like file (hg38 assets).
overlap_features_bed <- function(gr, bed_path, prefix) {
    p <- annot_path(bed_path)
    ann <- if (grepl("repeat-masker-hg38\\.gz$", p)) {
        dt <- fread(cmd = paste("zcat", shQuote(p)), header = FALSE, select = 1:3,
                    col.names = c("chrom", "start", "end"))
        GRanges(dt$chrom, IRanges(dt$start + 1L, dt$end))
    } else {
        rtracklayer::import(p)
    }
    overlap_features(gr, ann, prefix)
}

#' Broad genomic annotation: promoter / exonic / intronic / intergenic.
#'
#' Assigned by precedence, not by counting: a VMR overlapping both a promoter
#' and an intron is a promoter VMR. Precedence is fixed here rather than being
#' decided per run.
broad_genomic_annotation <- function(gr_hg38, gtf_bed, promoter_bp = 2000L) {
    dt <- fread(cmd = paste("cut -f1,2,3,6,8", shQuote(gtf_bed)), header = FALSE,
                col.names = c("chrom", "start", "end", "strand", "feature"))
    dt <- dt[feature %in% c("gene", "exon")]
    genes <- dt[feature == "gene"]
    exons <- dt[feature == "exon"]

    gene_gr <- GRanges(genes$chrom, IRanges(genes$start + 1L, genes$end),
                       strand = genes$strand)
    exon_gr <- GRanges(exons$chrom, IRanges(exons$start + 1L, exons$end))
    tss <- resize(gene_gr, width = 1L, fix = "start")
    prom_gr <- suppressWarnings(
        trim(promoters(tss, upstream = promoter_bp, downstream = promoter_bp)))

    ann <- rep("intergenic", length(gr_hg38))
    ann[overlapsAny(gr_hg38, gene_gr, ignore.strand = TRUE)] <- "intronic"
    ann[overlapsAny(gr_hg38, exon_gr, ignore.strand = TRUE)] <- "exonic"
    ann[overlapsAny(gr_hg38, prom_gr, ignore.strand = TRUE)] <- "promoter"
    factor(ann, levels = c("intergenic", "intronic", "exonic", "promoter"))
}
