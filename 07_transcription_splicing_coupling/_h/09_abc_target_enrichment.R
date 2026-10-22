#!/usr/bin/env Rscript
#### 07 stage 09 -- where do ABC target-gene sets fall along the gradient? ####
##
## Usage:
##   Rscript _h/09_abc_target_enrichment.R --cohort AA
##   Rscript _h/09_abc_target_enrichment.R --cohort AA --out-dir /path/to/draft
##   Rscript _h/09_abc_target_enrichment.R --cohort AA --allow-unlocked   # smoke
##
## Why this stage exists
## --------------------
## v1 connected intergenic VMRs to putative target genes through BrainScope
## Activity-by-Contact (ABC) enhancer-gene links and asked which pathways the
## target genes of heritable versus non-heritable VMRs fell into. v2 retired the
## heritable / non-heritable classes (AGENTS.md 7.2: local genetic control is
## continuous) and the pathway question was dropped with them. The PI amendment
## of 2026-10-22 (config/transcription_splicing.yml:abc_target_enrichment)
## restores it in rank form: each ABC target gene is scored by the mean
## standardized relative local genetic influence (z) of the VMRs sitting in its
## linked enhancers, and every gene set is asked whether its genes sit higher
## or lower on that gradient than the other ABC-linked genes.
##
## What it is NOT. It is not the expression_abc coupling modality (that tests
## VMR methylation against target-gene expression, and is excluded from the
## coupling FDR family on power grounds). No RNA-seq enters this stage. The
## universe is the set of genes the enhancers point at, whether or not they
## were quantified, because the question is about which genes the gradient
## indexes, not about coupling. Nothing here is a mediation, activity or
## target-gene regulation claim; the ABC links are a published resource, and
## a VMR inside a called enhancer is overlap, not function.
##
## Method
## ------
##   gene statistic   G_g = mean over linked VMRs v of z_v            (per cell)
##   set statistic    M_S = mean over genes g in S of G_g
##   null             permute z among the cell's scored VMRs, keeping every
##                    VMR -> gene link fixed; recompute G and M. Genes that share
##                    an enhancer VMR share a score, so a gene-level permutation
##                    would understate the null variance; a VMR-level one does
##                    not. Two-sided p with the +1 correction; BH within each
##                    gene-set collection within each cell.
##   pooled cell      the three regions' link sets concatenated, with z permuted
##                    within region so no region's scale leaks into another's.
##   secondary        fgsea multilevel on the same ranked gene list, reported
##                    beside the primary so a reader who expects a GSEA sees one;
##                    the share of tied statistics (shared VMRs) is recorded.
##
## Module-level: reads the accepted Module 07 runs and their accepted Module 02
## upstreams (one shared vmr_set_id per cell) and writes _m/combined/. It
## changes no decision. Exploratory, supplement-only (config:claims).

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

MODULE <- "07_transcription_splicing_coupling"
STAGE  <- "07_transcription_splicing_coupling/_h/09_abc_target_enrichment.R"
opts <- parse_v2_args(require = "cohort")
cohort <- opts$cohort
allow_unaccepted <- isTRUE(opts$allow_unlocked)

cfg_all <- load_config("transcription_splicing")
assert_locked(list(transcription_splicing = cfg_all), allow_unlocked = allow_unaccepted)
cfg <- cfg_all$abc_target_enrichment
if (!isTRUE(cfg$enabled)) stop("abc_target_enrichment is disabled in config")
stopifnot(identical(cfg$test, "vmr_permutation"),
          identical(cfg$gene_statistic, "mean_z"),
          identical(cfg$statistic, "mean_gene_z_in_set"))

regions <- as.character(unlist(load_config("cohorts")$regions))
cells <- as.character(unlist(cfg$cells))
stopifnot(all(setdiff(cells, "pooled") %in% regions))

module_root <- file.path(repo_root(), MODULE)
out_dir <- opts$out_dir %||% file.path(module_root, "_m", "combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sfx <- paste0("-", cohort)

n_perm <- as.integer(cfg$n_permutations)
seed   <- as.integer(cfg$seed)
alpha  <- as.numeric(cfg$fdr_threshold)
min_sz <- as.integer(cfg$min_set_size)
max_sz <- as.integer(cfg$max_set_size)

## ------------------------------------------------------------ upstreams
## The accepted Module 07 run names the cell; its accepted Module 02 upstream
## is the score of record. Both README rows carry vmr_set_id, and they must
## agree, or the links would be built on one catalog and scored on another.
tsc <- lapply(regions, function(r) require_accepted_upstream(
    MODULE, cohort, r, allow_unaccepted = allow_unaccepted))
lgv <- lapply(regions, function(r) require_accepted_upstream(
    "02_local_genetic_variance", cohort, r, allow_unaccepted = allow_unaccepted))
names(tsc) <- names(lgv) <- regions
smoke <- any(vapply(c(tsc, lgv), function(x) isTRUE(x$smoke_run), logical(1)))
for (r in regions) {
    if (!smoke && !identical(tsc[[r]]$vmr_set_id, lgv[[r]]$vmr_set_id)) {
        stop("vmr_set_id differs between accepted 07 and 02 runs for ", r, ": ",
             tsc[[r]]$vmr_set_id, " vs ", lgv[[r]]$vmr_set_id)
    }
}
if (smoke) warning("SMOKE: an upstream is unaccepted; outputs carry citable = FALSE")

scores <- rbindlist(lapply(regions, function(r) {
    d <- load_local_genetic_control(lgv[[r]]$run_id, r, cohort)
    d <- as.data.table(d)[local_genetic_control_eligible %in% c(TRUE, "TRUE")]
    d[, .(vmr_id, chrom, start, end, region = r,
          z = as.numeric(get(cfg$vmr_score)),
          score = as.numeric(local_snp_contribution_score))]
}))
if (anyNA(scores$z)) stop("NA in ", cfg$vmr_score, " among eligible VMRs")
message("[09] scored VMRs: ", paste(scores[, .N, by = region][, paste(region, N)], collapse = "; "))

## ------------------------------------------------------------ ABC links
## Same source, same overlap rule and same identifier handling as
## 01_build_feature_links.R's expression_abc branch, minus the assay-annotation
## merge (no expression filter here, see the header).
abc_f <- file.path(repo_root(), cfg$link_source)
if (!file.exists(abc_f)) stop("ABC link source not found: ", abc_f)
abc <- fread(abc_f)
need <- c("PeakID_seqnames", "PeakID_start", "PeakID_end", "TargetGene", "TargetGene_name")
if (!all(need %in% names(abc))) {
    stop("ABC file lacks expected columns: ", paste(setdiff(need, names(abc)), collapse = ", "))
}
norm_chr  <- function(x) { x <- as.character(x); ifelse(startsWith(x, "chr"), x, paste0("chr", x)) }
strip_ver <- function(x) sub("\\..*$", "", as.character(x))
abc2 <- unique(data.table(
    chrom = norm_chr(abc$PeakID_seqnames),
    start = as.integer(abc$PeakID_start),
    end   = as.integer(abc$PeakID_end),
    gene_id = strip_ver(abc$TargetGene),
    gene_name = as.character(abc$TargetGene_name)))
abc2 <- abc2[!is.na(gene_id) & nzchar(gene_id)]
egr <- GRanges(abc2$chrom, IRanges(abc2$start, abc2$end))

link_region <- function(r) {
    s <- scores[region == r]
    ## vmr.bed is 0-based half-open; vmr_id carries bed coordinates, and the
    ## range used for overlap is start + 1 .. end, as in 01_build_feature_links.R.
    vgr <- GRanges(s$chrom, IRanges(s$start + 1L, s$end))
    hits <- findOverlaps(vgr, egr)
    if (length(hits) == 0) stop("No VMR overlaps an ABC enhancer in ", r)
    unique(data.table(region = r,
                      vmr_id = s$vmr_id[queryHits(hits)],
                      z = s$z[queryHits(hits)],
                      gene_id = abc2$gene_id[subjectHits(hits)],
                      gene_name = abc2$gene_name[subjectHits(hits)]),
           by = c("region", "vmr_id", "gene_id"))
}
links <- rbindlist(lapply(regions, link_region))
message("[09] ABC links: ", paste(links[, .(n_vmr = uniqueN(vmr_id), n_gene = uniqueN(gene_id)),
                                        by = region][, paste0(region, " ", n_vmr, " VMRs / ", n_gene, " genes")],
                                  collapse = "; "))

## ------------------------------------------------------------ gene sets
suppressPackageStartupMessages(library(msigdbr))
## msigdbr renamed category/subcategory to collection/subcollection in 10.x;
## accept either so the stage runs against the Quest and local libraries.
msig_fetch <- function(collection, subcollection) {
    f <- names(formals(msigdbr))
    args <- if ("collection" %in% f) {
        list(species = "Homo sapiens", collection = collection, subcollection = subcollection)
    } else {
        list(species = "Homo sapiens", category = collection, subcategory = subcollection)
    }
    d <- as.data.table(do.call(msigdbr, args))
    ## 10.x exposes ensembl_gene for the target species; older releases name
    ## it the same way. Fail loudly rather than silently match zero genes.
    if (!"ensembl_gene" %in% names(d)) stop("msigdbr table lacks ensembl_gene")
    unique(d[!is.na(ensembl_gene) & nzchar(ensembl_gene),
             .(gs_name, gene_id = as.character(ensembl_gene))])
}
gs_tbl <- rbindlist(lapply(cfg$gene_sets, function(g) {
    d <- msig_fetch(g$collection, g$subcollection)
    d[, `:=`(collection = g$label, msig_collection = g$collection,
             msig_subcollection = g$subcollection)][]
}))
msig_version <- as.character(packageVersion("msigdbr"))
message("[09] gene sets: ", paste(gs_tbl[, uniqueN(gs_name), by = collection][, paste(collection, V1)],
                                   collapse = "; "), " (msigdbr ", msig_version, ")")

## ------------------------------------------------------------ the test
set.seed(seed)
`%by%` <- function(m, v) as.vector(m %*% v)

run_cell <- function(cell) {
    L <- if (cell == "pooled") copy(links) else links[region == cell]
    ## In the pooled cell a gene may be linked in several regions; it is one
    ## gene, scored over all its links, and z is permuted within region.
    L[, vkey := paste(region, vmr_id)]
    vmrs <- unique(L[, .(vkey, region, z)])
    genes <- unique(L[, .(gene_id, gene_name)])
    genes <- genes[!duplicated(gene_id)]
    n_v <- nrow(vmrs); n_g <- nrow(genes)
    vi <- match(L$vkey, vmrs$vkey); gi <- match(L$gene_id, genes$gene_id)
    ## Link matrix, genes x VMRs, rows normalised to the mean (gene_statistic).
    A <- matrix(0, n_g, n_v); A[cbind(gi, vi)] <- 1
    nl <- rowSums(A); stopifnot(all(nl >= as.integer(cfg$min_vmrs_per_gene)))
    A <- A / nl
    G_obs <- A %by% vmrs$z
    gene_scores <- data.table(cell = cell, gene_id = genes$gene_id,
                              gene_name = genes$gene_name, n_vmrs_linked = nl,
                              n_regions_linked = L[, uniqueN(region), by = gene_id][match(genes$gene_id, gene_id), V1],
                              mean_z = G_obs,
                              max_z = L[, max(z), by = gene_id][match(genes$gene_id, gene_id), V1])

    ## Gene sets restricted to the universe, then size-filtered within it.
    gs <- gs_tbl[gene_id %in% genes$gene_id]
    gs <- gs[, if (.N >= min_sz && .N <= max_sz) .SD, by = .(collection, gs_name)]
    if (nrow(gs) == 0) stop("No gene set within [", min_sz, ", ", max_sz, "] in ", cell)
    sets <- unique(gs[, .(collection, gs_name)])
    si <- match(paste(gs$collection, gs$gs_name), paste(sets$collection, sets$gs_name))
    S <- matrix(0, nrow(sets), n_g); S[cbind(si, match(gs$gene_id, genes$gene_id))] <- 1
    sz <- rowSums(S); S <- S / sz
    M_obs <- S %by% G_obs

    ## Null: permute z among VMRs within region (a single region is one block).
    blocks <- split(seq_len(n_v), vmrs$region)
    Zp <- matrix(NA_real_, n_v, n_perm)
    for (b in blocks) Zp[b, ] <- replicate(n_perm, vmrs$z[b][sample.int(length(b))])
    M_perm <- S %*% (A %*% Zp)                       # sets x n_perm
    mu <- rowMeans(M_perm)
    p_perm <- (rowSums(abs(M_perm - mu) >= abs(M_obs - mu)) + 1) / (n_perm + 1)
    res <- data.table(cell = cell, collection = sets$collection, gs_name = sets$gs_name,
                      n_genes_in_set = as.integer(sz), n_genes_universe = n_g,
                      n_vmrs_universe = n_v, mean_z_in_set = M_obs,
                      mean_z_universe = mean(G_obs), delta_z = M_obs - mean(G_obs),
                      perm_null_mean = mu, perm_null_sd = apply(M_perm, 1, sd),
                      p_perm = p_perm, n_permutations = n_perm)
    res[, q_bh := p.adjust(p_perm, method = "BH"), by = collection]
    res[, significant := q_bh < alpha]

    ## Secondary: fgsea on the same ranked list. Ties come from genes sharing
    ## a VMR; their share is recorded so the warning is a number, not a surprise.
    if (isTRUE(cfg$secondary_fgsea)) {
        suppressPackageStartupMessages(library(fgsea))
        stat <- setNames(G_obs, genes$gene_id)
        tie_frac <- 1 - length(unique(stat)) / length(stat)
        pathways <- split(gs$gene_id, paste(gs$collection, gs$gs_name, sep = "||"))
        fg <- suppressWarnings(as.data.table(fgseaMultilevel(
            pathways, stat, minSize = min_sz, maxSize = max_sz, eps = 0,
            scoreType = "std", nPermSimple = 2000)))
        fg[, c("collection", "gs_name") := tstrsplit(pathway, "||", fixed = TRUE)]
        fg[, fgsea_q_bh := p.adjust(pval, method = "BH"), by = collection]
        res <- merge(res, fg[, .(collection, gs_name, fgsea_nes = NES, fgsea_p = pval,
                                 fgsea_q_bh, fgsea_leading_edge_n = lengths(leadingEdge))],
                     by = c("collection", "gs_name"), all.x = TRUE)
        res[, fgsea_tie_fraction := tie_frac]
    }
    setorder(res, collection, p_perm)
    list(res = res, genes = gene_scores)
}

out <- lapply(cells, run_cell)
enr   <- rbindlist(lapply(out, `[[`, "res"), fill = TRUE)
genes <- rbindlist(lapply(out, `[[`, "genes"))

## ------------------------------------------------------------ outputs
citable <- !smoke
tag <- function(d) {
    d[, `:=`(cohort = cohort, citable = citable,
             exploratory_supplement_only = isTRUE(cfg$exploratory_supplement_only),
             post_hoc_addition = isTRUE(cfg$post_hoc_addition),
             causal_claim_allowed = FALSE, stage = STAGE)][]
}
enr   <- tag(enr)
genes <- tag(genes)
links_out <- tag(links[, .(region, vmr_id, z, gene_id, gene_name)])

decision <- enr[, .(n_sets_tested = .N, n_sets_q_lt_alpha = sum(significant),
                    n_sets_higher = sum(significant & delta_z > 0),
                    n_sets_lower  = sum(significant & delta_z < 0),
                    min_q = min(q_bh), n_genes_universe = n_genes_universe[1],
                    n_vmrs_universe = n_vmrs_universe[1]),
                by = .(cell, collection)]
decision[, reading := fifelse(n_sets_q_lt_alpha == 0,
                              "NO_SET_SEPARATED_ALONG_GRADIENT",
                              "SETS_SEPARATED_ALONG_GRADIENT_EXPLORATORY")]
decision <- tag(decision)

prov <- data.table(
    field = c("stage", "cohort", "run_date", "source_07_runs", "source_02_runs",
              "vmr_set_ids", "abc_link_source", "abc_link_source_sha256",
              "config_sha256", "msigdbr_version", "msigdb_db_version",
              "n_permutations", "seed", "gene_sets", "min_set_size", "max_set_size",
              "fdr", "git_commit", "git_dirty", "r_version", "citable"),
    value = c(STAGE, cohort, format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
              paste(vapply(tsc, `[[`, "", "run_id"), collapse = ";"),
              paste(vapply(lgv, `[[`, "", "run_id"), collapse = ";"),
              paste(vapply(tsc, function(x) x$vmr_set_id %||% NA_character_, ""), collapse = ";"),
              cfg$link_source, file_sha256(abc_f),
              attr(cfg_all, "config_sha256"), msig_version,
              tryCatch(as.character(unique(as.data.table(msigdbr(species = "Homo sapiens",
                       collection = "H"))$db_version))[1], error = function(e) NA_character_),
              n_perm, seed,
              paste(vapply(cfg$gene_sets, function(g) paste0(g$collection, "/", g$subcollection), ""),
                    collapse = ";"),
              min_sz, max_sz, paste(cfg$fdr_method, cfg$fdr_family, alpha),
              tryCatch(system2(Sys.getenv("V2_GIT_BIN", "git"),
                               c("-C", shQuote(repo_root()), "rev-parse", "HEAD"),
                               stdout = TRUE), error = function(e) NA_character_),
              tryCatch(as.character(git_dirty(repo_root())), error = function(e) NA_character_),
              R.version.string, citable))

write_atomic(enr,       file.path(out_dir, paste0("abc-target-enrichment", sfx, ".tsv")))
write_atomic(genes,     file.path(out_dir, paste0("abc-target-gene-scores", sfx, ".tsv")))
write_atomic(links_out, file.path(out_dir, paste0("abc-target-links", sfx, ".tsv")))
write_atomic(decision,  file.path(out_dir, paste0("abc-target-enrichment-decision", sfx, ".tsv")))
write_atomic(prov,      file.path(out_dir, paste0("abc-target-enrichment-provenance", sfx, ".tsv")))

message("[09] done. ", paste(decision[, paste0(cell, "/", collection, ": ",
                                              n_sets_q_lt_alpha, "/", n_sets_tested,
                                              " sets q<", alpha)], collapse = "; "))
