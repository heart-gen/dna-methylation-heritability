#!/usr/bin/env Rscript
#### 09_schizophrenia_gwas_loci -- EA meQTL x PGC3 European coloc ####
##
## Usage (coloc env, 00_shared/slurm.sh:run_r_coloc):
##   Rscript _h/21_ea_coloc.R --run-id scz-all_individuals.EA-crossregion-YYYYMMDD
##
## EXPLORATORY SUPPORT, NOT A GATE (config/scz_ea_targeted.yml:coloc.gating).
## Eligibility is fixed in config before any EA result existed: the locus's
## lead pair reproduces in EA, and the lead CpG's EA cis scan reaches
## eligibility_min_p; at most max_loci, strongest first, across regions.
##
## This is the LD-matched pairing the AA meQTL arm (Module 09 arm_b) could never
## be: PGC3 European GWAS x a meQTL scan in white-American donors. coloc.abf's
## single-causal-variant model still applies, and the EA cells are small
## (55-129 donors), so a low PP4 here is weak evidence against sharing.

suppressPackageStartupMessages({
    library(data.table)
    library(coloc)
    library(yaml)
})
fread_gz <- function(f, ...) {
    if (grepl("\\.gz$", f)) data.table::fread(cmd = paste("zcat", shQuote(f)), ...)
    else data.table::fread(f, ...)
}
repo_root <- function(start = getwd()) {
    dir <- normalizePath(Sys.getenv("V2_REPO_ROOT", start), mustWork = TRUE)
    while (dir != dirname(dir)) {
        if (dir.exists(file.path(dir, ".git"))) return(dir)
        dir <- dirname(dir)
    }
    stop("Could not locate repository root")
}
args <- commandArgs(trailingOnly = TRUE)
run_id <- args[which(args == "--run-id") + 1]
MODULE <- "09_schizophrenia_gwas_loci"
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", run_id)
cfg <- read_yaml(file.path(run_dir, "code", "config", "scz_ea_targeted.yml"))
scz <- read_yaml(file.path(run_dir, "code", "config", "schizophrenia.yml"))
cc <- cfg$coloc
regions <- names(cfg$accepted_runs)

scans <- rbindlist(lapply(regions, function(re) {
    f <- file.path(run_dir, "results", re, "cis-scan-summary.tsv")
    if (!file.exists(f) || file.size(f) == 0) return(NULL)
    d <- fread(f)
    if (!nrow(d)) return(NULL)
    d[, region := re][]
}), fill = TRUE)

out_f <- file.path(run_dir, "results", "ea-coloc.tsv")
if (!nrow(scans)) {
    fwrite(data.table(status = "NO_REPRODUCING_LOCUS",
                      message = "no lead pair reproduced in EA; nothing is coloc-eligible",
                      gating = FALSE), out_f, sep = "\t")
    quit(save = "no")
}
scans[, eligible := min_pval <= cc$eligibility_min_p]
setorder(scans, -eligible, min_pval)
scans[, selected := eligible & seq_len(.N) <= cc$max_loci]

n_gwas <- scz$locus_definition$gwas_n_cases + scz$locus_definition$gwas_n_controls
s_gwas <- scz$locus_definition$gwas_n_cases / n_gwas

res <- rbindlist(lapply(seq_len(nrow(scans)), function(i) {
    r <- scans[i]
    base <- data.table(region = r$region, locus_id = r$locus_id, cpg_id = r$lead_cpg_id,
                       ea_cis_min_p = r$min_pval, eligible = r$eligible, selected = r$selected,
                       n_variants_shared = NA_integer_, PP0 = NA_real_, PP1 = NA_real_,
                       PP2 = NA_real_, PP3 = NA_real_, PP4 = NA_real_,
                       best_gwas_pvalue = NA_real_, status = NA_character_, gating = FALSE)
    if (!r$selected) {
        base[, status := if (r$eligible) "ELIGIBLE_NOT_SELECTED_MAX_LOCI" else "NOT_ELIGIBLE_WEAK_EA_SIGNAL"]
        return(base)
    }
    q <- fread_gz(file.path(run_dir, "results", r$region, paste0("cis-scan-locus", r$locus_id, ".tsv.gz")))
    ## other allele = the ID allele that is not counted
    q[, other_allele := fifelse(counted_allele == id_allele_1, id_allele_2, id_allele_1)]
    q[, key := paste0(chrom, "_", pos)]
    g <- fread_gz(file.path(repo_root(), MODULE, "_m", "runs", cfg$accepted_runs[[r$region]],
                            "gwas", "gwas-sumstats-hg38.tsv.gz"))
    cpos <- as.integer(sub(".*:", "", r$lead_cpg_id))
    g <- g[chrom == q$chrom[1] & abs(pos_hg38 - cpos) <= cc$window_bp]
    g[, key := paste0(chrom, "_", pos_hg38)]
    m <- merge(q, g, by = "key")
    m[, same := effect_allele == counted_allele & other_allele.y == other_allele.x]
    m[, flip := effect_allele == other_allele.x & other_allele.y == counted_allele]
    m <- m[same | flip]
    m[, gwas_beta := fifelse(same, beta, -beta)]
    m <- m[is.finite(slope) & is.finite(slope_se) & slope_se > 0 & is.finite(se) & se > 0]
    m[, maf := pmin(counted_freq, 1 - counted_freq)]
    m <- unique(m[maf > 0 & maf < 0.5], by = "key")
    base[, n_variants_shared := nrow(m)]
    if (nrow(m) < cc$min_variants_shared) {
        base[, status := "UNEVALUABLE_TOO_FEW_SHARED_VARIANTS"]
        return(base)
    }
    d1 <- list(beta = m$gwas_beta, varbeta = m$se^2, snp = m$key, type = "cc",
               s = s_gwas, N = n_gwas, MAF = m$maf)
    d2 <- list(beta = m$slope, varbeta = m$slope_se^2, snp = m$key, type = "quant",
               N = max(m$n), MAF = m$maf)
    fit <- suppressWarnings(coloc.abf(d1, d2, p1 = cc$p1, p2 = cc$p2, p12 = cc$p12))
    s <- fit$summary
    base[, `:=`(PP0 = s[["PP.H0.abf"]], PP1 = s[["PP.H1.abf"]], PP2 = s[["PP.H2.abf"]],
                PP3 = s[["PP.H3.abf"]], PP4 = s[["PP.H4.abf"]],
                best_gwas_pvalue = min(m$pvalue), status = "OK")]
    base
}), fill = TRUE)
res[, pp4_meets_threshold := !is.na(PP4) & PP4 >= cc$pp4_threshold]
fwrite(res, out_f, sep = "\t")
print(res)
