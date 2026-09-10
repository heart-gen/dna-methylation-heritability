#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- coloc.susie sensitivity ####
##
## Usage:
##   Rscript _h/11_coloc_susie_sensitivity.R --run-id scz-AA-caudate-YYYYMMDD
##
## Runs in the coloc conda environment (00_shared/slurm.sh:run_r_coloc).
##
## coloc.abf assumes AT MOST ONE causal variant per region. Schizophrenia loci
## frequently carry several, and where they do, a single-causal-variant model
## can report H3 (distinct signals) for a region in which one of several shared
## signals is genuinely shared. This stage is the prespecified sensitivity check
## on that assumption, restricted to the prioritized loci because it is the
## prioritized loci that reach the manuscript.
##
## It uses the credible sets GTEx v11 already ships
## (*.eQTLs.SuSiE_summary.parquet) rather than refitting SuSiE on the QTL side:
## those were fitted on the individual-level genotypes with the true LD matrix,
## which is strictly better than anything reconstructible from summary
## statistics here. What that gives is a credible-set-level answer -- does a
## GTEx credible set for this gene contain the variants that carry the GWAS
## signal -- and it is reported as such.
##
## This is a SENSITIVITY analysis. It never overturns the coloc.abf result and
## never sets claimable_colocalization; it is read alongside it. A locus where
## the two disagree is reported as disagreeing.

suppressPackageStartupMessages({
    library(data.table)
    library(arrow)
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

root <- repo_root()
run_dir <- file.path(root, "09_schizophrenia_risk_application", "_m", "runs", run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

cfg_f <- file.path(run_dir, "code", "config", "schizophrenia.yml")
if (!file.exists(cfg_f)) cfg_f <- file.path(root, "config", "schizophrenia.yml")
scz <- yaml::read_yaml(cfg_f)
gtex <- scz$colocalization$gtex

out_f <- file.path(run_dir, "results", "coloc-susie.tsv")
prio_f <- file.path(run_dir, "results", "prioritized-loci.tsv")
if (!file.exists(prio_f)) stop("Run 10_prioritize_loci.R first: ", prio_f)
prio <- fread(prio_f)
prio <- if (nrow(prio)) prio[prioritized == TRUE] else prio

if (nrow(prio) == 0) {
    fwrite(data.table(), out_f, sep = "\t")
    message("[09] no prioritized locus; SuSiE sensitivity not run")
    quit(status = 0)
}

abf_f <- file.path(run_dir, "results", "coloc-abf.tsv.gz")
if (!file.exists(abf_f)) stop("Run 09b_combine_coloc.R first: ", abf_f)
abf <- fread_gz(abf_f)
## Sensitivity applies only to the gate-eligible, ancestry-matched arms; the
## cross-ancestry arm has a prior problem that a multi-causal-variant model
## does not fix.
abf <- abf[locus_id %in% prio$locus_id & gate_eligible == TRUE &
           ld_ancestry_matched == TRUE & status == "OK"]
if (nrow(abf) == 0) {
    fwrite(data.table(), out_f, sep = "\t")
    message("[09] no gate-eligible evaluated region among prioritized loci")
    quit(status = 0)
}

region_dir <- file.path(run_dir, "coloc", "regions")
variant_key <- function(v) {
    p <- tstrsplit(as.character(v), "_", fixed = TRUE)
    a <- toupper(p[[3]]); b <- toupper(p[[4]])
    paste0(p[[1]], "_", p[[2]], "_", pmin(a, b), "_", pmax(a, b))
}

susie_file <- function(arm, tissue) {
    dir <- if (identical(arm, "gtex_sqtl")) gtex$sqtl_dir else gtex$eqtl_dir
    kind <- if (identical(arm, "gtex_sqtl")) "sQTLs" else "eQTLs"
    file.path(dir, paste0(tissue, ".v11.", kind, ".SuSiE_summary.parquet"))
}

cs_cache <- new.env(parent = emptyenv())
load_cs <- function(arm, tissue) {
    key <- paste(arm, tissue, sep = "|")
    if (!is.null(cs_cache[[key]])) return(cs_cache[[key]])
    f <- susie_file(arm, tissue)
    cs <- if (file.exists(f)) {
        as.data.table(arrow::read_parquet(
            f, col_select = c("phenotype_id", "variant_id", "pip", "cs_id",
                              "cs_size")))
    } else data.table()
    cs_cache[[key]] <- cs
    cs
}

rows <- lapply(seq_len(nrow(abf)), function(i) {
    r <- as.list(abf[i])
    base <- data.table(
        run_id = run_id, locus_id = r$locus_id, arm = r$arm,
        qtl_context = r$qtl_context, phenotype_id = r$phenotype_id,
        abf_pp4 = r$PP4, abf_colocalized = as.logical(r$colocalized),
        n_credible_sets = NA_integer_, cs_id = NA_integer_, cs_size = NA_integer_,
        n_cs_variants_in_region = NA_integer_,
        max_pip_among_gwas_lead = NA_real_,
        gwas_lead_in_credible_set = NA, status = "NOT_RUN", message = "")

    cs <- load_cs(r$arm, r$qtl_context)
    if (nrow(cs) == 0) {
        base[, `:=`(status = "NO_SUSIE_SUMMARY",
                    message = susie_file(r$arm, r$qtl_context))]
        return(base)
    }
    cs <- cs[phenotype_id == r$phenotype_id]
    if (nrow(cs) == 0) {
        base[, `:=`(status = "NO_CREDIBLE_SET_FOR_PHENOTYPE")]
        return(base)
    }
    f <- file.path(region_dir, paste0(
        r$chrom, "__", r$locus_id, "__", r$arm, "__", r$qtl_context, "__",
        gsub("/", "_", r$phenotype_id), ".tsv.gz"))
    if (!file.exists(f)) {
        base[, `:=`(status = "REGION_FILE_MISSING", message = f)]
        return(base)
    }
    d <- fread_gz(f)
    cs[, match_key := variant_key(variant_id)]
    cs <- cs[match_key %in% d$match_key]
    if (nrow(cs) == 0) {
        base[, `:=`(status = "CREDIBLE_SET_OUTSIDE_REGION",
                    n_credible_sets = 0L)]
        return(base)
    }
    ## "Does the QTL credible set contain the variants carrying the GWAS
    ## signal?" The GWAS lead is the strongest GWAS variant of the harmonised
    ## region, so this is a statement about shared variants, not about which
    ## variant is causal.
    lead <- d[which.min(gwas_pvalue)]
    best <- cs[, .(max_pip = max(pip, na.rm = TRUE),
                   cs_size = max(cs_size, na.rm = TRUE),
                   n_in_region = .N,
                   contains_lead = any(match_key == lead$match_key)),
               by = cs_id]
    setorder(best, -contains_lead, -max_pip)
    top <- best[1]
    base[, `:=`(n_credible_sets = nrow(best), cs_id = top$cs_id,
                cs_size = top$cs_size, n_cs_variants_in_region = top$n_in_region,
                max_pip_among_gwas_lead =
                    cs[match_key == lead$match_key, if (.N) max(pip) else NA_real_],
                gwas_lead_in_credible_set = top$contains_lead,
                status = "OK")]
    base
})

out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
out[, agrees_with_abf := status == "OK" &
        (as.logical(gwas_lead_in_credible_set) == abf_colocalized)]
out[, method := "coloc.susie_credible_set_overlap"]
out[, role := "sensitivity_only_never_sets_a_claim"]
setorder(out, locus_id, arm, qtl_context, phenotype_id)
# write_atomic() lives in 00_shared/load.R, which this stage does not source:
# it runs in the coloc env, whose R has neither here() nor the shared helpers.
# Same tmp-then-rename, inlined, so a killed job cannot leave a partial table.
tmp_f <- paste0(out_f, ".tmp")
fwrite(out, tmp_f, sep = "\t")
if (!file.rename(tmp_f, out_f)) stop("could not rename ", tmp_f, " -> ", out_f)

message("[09] SuSiE sensitivity: ", sum(out$status == "OK"), "/", nrow(out),
        " regions evaluated; ", sum(out$agrees_with_abf, na.rm = TRUE),
        " agree with coloc.abf")
