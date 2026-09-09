#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- prioritize illustrative loci ####
##
## Usage:
##   Rscript _h/10_prioritize_loci.R --run-id scz-AA-caudate-YYYYMMDD
##
## AGENTS.md 7.8: "prioritize at most five illustrative loci using prespecified
## criteria", and "existing hero loci rs8048039 and rs13331198 remain
## candidates, not fixed answers. Reprioritize after the corrected run." The
## rule is config/schizophrenia.yml:prioritization.rule (ranked_composite_v1),
## locked by the PI on 2026-09-09 BEFORE this run existed. This script applies
## it mechanically; it contains no locus names.
##
## Rule, in order:
##   filter   the locus has >= 1 FDR-significant risk-variant x CpG meQTL link
##   rank 1   max |local_snp_contribution_score_z| over its linked VMRs, desc
##   rank 2   has transcriptional coupling (Module 07), desc
##   rank 3   PGC3 index-SNP p-value, asc
##   tie      locus_id, asc  -- so the output is byte-stable across reruns
##
## Every ranking key is written out beside the rank, so the ordering can be
## re-derived by hand from this one table. Colocalization is deliberately NOT a
## ranking key: it is downstream evidence about the prioritized loci, and using
## it to choose them would make the coloc result partly self-selected.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09_schizophrenia_risk_application", "_h")),
    "run_config.R"))

suppressPackageStartupMessages(library(data.table))

MODULE <- "09_schizophrenia_risk_application"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mf("cohort"); region <- mf("region")
scz <- load_run_config("schizophrenia", run_dir)
rule <- scz$prioritization$rule
if (!identical(rule$name, "ranked_composite_v1")) {
    stop("Unrecognised prioritization rule '", rule$name, "'")
}
max_loci <- as.integer(scz$prioritization$max_loci)

read_or_empty <- function(name) {
    f <- file.path(run_dir, "results", name)
    if (file.exists(f)) fread(f) else data.table()
}

loci <- fread(file.path(run_dir, "results", "scz-loci.tsv"))
support <- read_or_empty("locus-meqtl-support.tsv")
coupling <- read_or_empty("locus-coupling.tsv")
gtex <- read_or_empty("gtex-support.tsv")
coloc <- read_or_empty("coloc-locus-gate-eligible.tsv")
links <- fread(file.path(run_dir, "results", "locus-vmr-links.tsv"))

## Ranking key 1, computed here from the linked VMRs so it is reproducible from
## a single table rather than from a filter applied in an earlier script.
axis <- links[, .(max_abs_local_snp_contribution_score_z =
                      max(abs(local_snp_contribution_score_z), na.rm = TRUE),
                  n_linked_vmrs = uniqueN(vmr_id)),
              by = locus_id]

dt <- loci[, .(locus_id, chrom, start, end, index_snp, pgc3_index_pvalue)]
dt <- merge(dt, axis, by = "locus_id", all.x = TRUE)
if (nrow(support)) {
    dt <- merge(dt, support[, .(locus_id, n_pairs_tested, n_pairs_significant,
                                n_vmrs_significant, min_qvalue,
                                has_significant_risk_variant_cpg_meqtl)],
                by = "locus_id", all.x = TRUE)
} else {
    dt[, `:=`(n_pairs_tested = 0L, n_pairs_significant = 0L,
              n_vmrs_significant = 0L, min_qvalue = NA_real_,
              has_significant_risk_variant_cpg_meqtl = FALSE)]
}
if (nrow(coupling)) {
    dt <- merge(dt, coupling[, .(locus_id, n_vmrs_coupled_any_modality,
                                 modalities_coupled,
                                 has_transcriptional_coupling)],
                by = "locus_id", all.x = TRUE)
} else {
    dt[, `:=`(n_vmrs_coupled_any_modality = 0L, modalities_coupled = "",
              has_transcriptional_coupling = FALSE)]
}
if (nrow(gtex)) {
    dt <- merge(dt, gtex[, .(locus_id, n_gtex_genes, gtex_tissues,
                             has_external_genetic_support)],
                by = "locus_id", all.x = TRUE)
} else {
    dt[, `:=`(n_gtex_genes = 0L, gtex_tissues = "",
              has_external_genetic_support = FALSE)]
}
if (nrow(coloc)) {
    dt <- merge(dt, coloc, by = "locus_id", all.x = TRUE)
} else {
    dt[, `:=`(max_pp4_gate_eligible = NA_real_,
              has_claimable_colocalization = FALSE)]
}

for (c_ in c("has_significant_risk_variant_cpg_meqtl",
             "has_transcriptional_coupling", "has_external_genetic_support",
             "has_claimable_colocalization")) {
    dt[[c_]] <- !is.na(dt[[c_]]) & as.logical(dt[[c_]])
}
for (c_ in c("n_pairs_tested", "n_pairs_significant", "n_vmrs_significant",
             "n_vmrs_coupled_any_modality", "n_gtex_genes", "n_linked_vmrs")) {
    dt[is.na(get(c_)), (c_) := 0L]
}

## ------------------------------------------------------------------ the rule
eligible <- dt[has_significant_risk_variant_cpg_meqtl == TRUE]
## A missing score cannot outrank a real one; -Inf sorts it last under desc
## rather than letting NA float to the top.
eligible[, rank_key_1 := fifelse(is.finite(max_abs_local_snp_contribution_score_z),
                                 max_abs_local_snp_contribution_score_z, -Inf)]
eligible[, rank_key_2 := as.integer(has_transcriptional_coupling)]
eligible[, rank_key_3 := fifelse(is.finite(pgc3_index_pvalue),
                                 pgc3_index_pvalue, Inf)]
setorder(eligible, -rank_key_1, -rank_key_2, rank_key_3, locus_id)
eligible[, rank := seq_len(.N)]
eligible[, prioritized := rank <= max_loci]

eligible[, `:=`(run_id = opts$run_id, cohort = cohort, region = region,
                rule_name = rule$name, max_loci = max_loci)]
setcolorder(eligible, c("rank", "prioritized", "locus_id", "chrom", "start",
                        "end", "index_snp", "rank_key_1", "rank_key_2",
                        "rank_key_3"))
write_atomic(eligible, file.path(run_dir, "results", "prioritized-loci.tsv"))

## The full evidence table for every locus, prioritized or not, so the five are
## readable in context and a reviewer can check what was left out.
## copy() first: `dt` is the result of a chain of merges, so := on it warns
## about taking a shallow copy.
dt <- copy(dt)
dt[, prioritized := locus_id %in% eligible[prioritized == TRUE, locus_id]]
setorder(dt, -prioritized, -has_significant_risk_variant_cpg_meqtl,
         min_qvalue, locus_id)
write_atomic(dt, file.path(run_dir, "results", "locus-evidence.tsv"))

message("[09] ", nrow(eligible), " eligible loci; prioritized ",
        sum(eligible$prioritized), " (cap ", max_loci, ")")
if (nrow(eligible)) {
    print(eligible[prioritized == TRUE,
                   .(rank, locus_id, index_snp, rank_key_1, rank_key_2, rank_key_3)])
}
