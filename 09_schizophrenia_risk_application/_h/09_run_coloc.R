#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- coloc.abf over prepared regions ####
##
## Usage (inside an array task):
##   Rscript _h/09_run_coloc.R --run-id scz-AA-caudate-YYYYMMDD --chrom 22
##
## Runs in the coloc conda environment (00_shared/slurm.sh:run_r_coloc); the
## epigenomics env has neither coloc nor arrow.
##
## AGENTS.md 7.8 permits a colocalization CLAIM only when the analysis "has been
## run with adequate ancestry-matched LD and passes its own gate". That
## condition is a property of the ARM, not of the posterior, so it is decided
## here per region and carried in the output rather than argued afterwards:
##
##   arm_a (GTEx eQTL/sQTL)  PGC3 European x GTEx European. LD matched.
##                           Gate-eligible; may support a claim.
##   arm_b (CpG meQTL)       PGC3 European GWAS x African-American meQTL.
##                           coloc.abf assumes both studies share an LD
##                           structure, and they do not. Recorded, never
##                           claimable, and excluded from the gate.
##
## Two other things are deliberate. Priors come from config, so a rerun cannot
## drift onto different defaults. And a region with too few shared variants is
## marked UNEVALUABLE rather than assigned a posterior: coloc.abf will happily
## return PP4 = 0.99 from twenty variants, and that number would be indefensible.

suppressPackageStartupMessages({
    library(data.table)
    library(coloc)
    library(yaml)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

## data.table::fread() reads .gz only when R.utils is installed, and it is not
## in the coloc environment. Piping through zcat keeps this stage working
## without adding a package to a shared env that other analyses also use.
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
get_arg <- function(flag, default = NULL) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) return(default)
    args[[i + 1]]
}
run_id <- get_arg("--run-id") %||% stop("--run-id is required")
chrom_arg <- get_arg("--chrom") %||% stop("--chrom is required")
chrom <- paste0("chr", sub("^chr", "", chrom_arg))

root <- repo_root()
run_dir <- file.path(root, "09_schizophrenia_risk_application", "_m", "runs", run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

cfg_f <- file.path(run_dir, "code", "config", "schizophrenia.yml")
if (!file.exists(cfg_f)) cfg_f <- file.path(root, "config", "schizophrenia.yml")
scz <- yaml::read_yaml(cfg_f)
cc <- scz$colocalization

region_dir <- file.path(run_dir, "coloc", "regions")
index_f <- file.path(region_dir, paste0(chrom, ".index.tsv"))
if (!file.exists(index_f)) {
    stop("No prepared-region index for ", chrom, "; run 08 first: ", index_f)
}
idx <- fread(index_f)

out_dir <- file.path(run_dir, "coloc", "abf")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_f <- file.path(out_dir, paste0(chrom, ".tsv.gz"))

if (nrow(idx) == 0) {
    fwrite(data.table(), out_f, sep = "\t", compress = "gzip")
    message("[09] ", chrom, ": no prepared regions")
    quit(status = 0)
}

min_shared <- as.numeric(cc$min_variants_shared)
p1 <- as.numeric(cc$p1); p2 <- as.numeric(cc$p2); p12 <- as.numeric(cc$p12)

## coloc's own input checks are strict and informative; a region that fails one
## is recorded with its message rather than aborting the chromosome. Losing 200
## regions because the 3rd was malformed would be the worse failure.
run_one <- function(r) {
    base <- data.table(
        chrom = r$chrom, locus_id = r$locus_id, arm = r$arm,
        qtl_context = r$qtl_context, phenotype_id = r$phenotype_id,
        n_variants_shared = r$n_variants_shared,
        qtl_ancestry = r$qtl_ancestry,
        ld_ancestry_matched = as.logical(r$ld_ancestry_matched),
        gate_eligible = as.logical(r$gate_eligible),
        arm_status = r$arm_status,
        PP0 = NA_real_, PP1 = NA_real_, PP2 = NA_real_, PP3 = NA_real_,
        PP4 = NA_real_, pp4_over_pp3 = NA_real_,
        best_gwas_pvalue = NA_real_, best_qtl_pvalue = NA_real_,
        status = NA_character_, message = "")

    if (!isTRUE(r$n_variants_shared >= min_shared) || !nzchar(r$region_file)) {
        base[, `:=`(status = "UNEVALUABLE_TOO_FEW_SHARED_VARIANTS",
                    message = paste0(r$n_variants_shared, " shared variants < ",
                                     min_shared))]
        return(base)
    }
    f <- file.path(region_dir, r$region_file)
    if (!file.exists(f)) {
        base[, `:=`(status = "UNEVALUABLE_REGION_FILE_MISSING", message = f)]
        return(base)
    }
    d <- fread_gz(f)
    d <- d[is.finite(gwas_beta) & is.finite(gwas_se) & gwas_se > 0 &
           is.finite(qtl_beta) & is.finite(qtl_se) & qtl_se > 0]
    ## MAF must be strictly inside (0, 0.5]; coloc rejects 0 and 1, and an
    ## allele frequency reported as > 0.5 is the frequency of the other allele.
    d[, maf := pmin(qtl_af, 1 - qtl_af)]
    d <- d[is.finite(maf) & maf > 0 & maf < 0.5]
    d <- unique(d, by = "match_key")
    if (nrow(d) < min_shared) {
        base[, `:=`(status = "UNEVALUABLE_TOO_FEW_SHARED_VARIANTS",
                    n_variants_shared = nrow(d),
                    message = "after finite-statistic and MAF filtering")]
        return(base)
    }

    ds_gwas <- list(beta = d$gwas_beta, varbeta = d$gwas_se^2,
                    snp = d$match_key, type = r$gwas_type,
                    N = as.numeric(r$gwas_n), MAF = d$maf)
    if (identical(r$gwas_type, "cc")) ds_gwas$s <- as.numeric(r$gwas_s)
    ds_qtl <- list(beta = d$qtl_beta, varbeta = d$qtl_se^2,
                   snp = d$match_key, type = "quant",
                   N = as.numeric(r$qtl_n), MAF = d$maf)

    res <- tryCatch(
        suppressWarnings(coloc::coloc.abf(ds_gwas, ds_qtl,
                                          p1 = p1, p2 = p2, p12 = p12)),
        error = function(e) e)
    if (inherits(res, "error")) {
        base[, `:=`(status = "FAILED", message = conditionMessage(res),
                    n_variants_shared = nrow(d))]
        return(base)
    }
    s <- res$summary
    base[, `:=`(
        n_variants_shared = nrow(d),
        PP0 = unname(s["PP.H0.abf"]), PP1 = unname(s["PP.H1.abf"]),
        PP2 = unname(s["PP.H2.abf"]), PP3 = unname(s["PP.H3.abf"]),
        PP4 = unname(s["PP.H4.abf"]),
        pp4_over_pp3 = unname(s["PP.H4.abf"]) / unname(s["PP.H3.abf"]),
        best_gwas_pvalue = min(d$gwas_pvalue, na.rm = TRUE),
        best_qtl_pvalue = min(d$qtl_pvalue, na.rm = TRUE),
        status = "OK")]
    base
}

results <- rbindlist(lapply(seq_len(nrow(idx)), function(i) run_one(as.list(idx[i]))),
                     use.names = TRUE, fill = TRUE)

thr <- as.numeric(cc$pp4_threshold)
ratio_min <- as.numeric(cc$pp4_over_pp3_min)
results[, colocalized := status == "OK" & is.finite(PP4) & PP4 >= thr &
                         is.finite(pp4_over_pp3) & pp4_over_pp3 >= ratio_min]
## The single most important column in this file. A cross-ancestry arm can
## produce a high PP4 and it still may not be called colocalization.
results[, claimable_colocalization := colocalized & gate_eligible &
                                      ld_ancestry_matched]
results[, `:=`(run_id = run_id, method = "coloc.abf",
               pp4_threshold = thr, pp4_over_pp3_min = ratio_min,
               p1 = p1, p2 = p2, p12 = p12)]
setorder(results, -PP4, locus_id, arm)
fwrite(results, out_f, sep = "\t", compress = "gzip")

message("[09] ", chrom, ": ", nrow(results), " regions, ",
        sum(results$status == "OK"), " evaluated, ",
        sum(results$claimable_colocalization, na.rm = TRUE), " claimable")
