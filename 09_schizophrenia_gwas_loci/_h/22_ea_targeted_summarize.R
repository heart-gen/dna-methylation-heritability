#!/usr/bin/env Rscript
#### 09_schizophrenia_gwas_loci -- targeted EA check: read and seal ####
##
## Usage:
##   Rscript _h/22_ea_targeted_summarize.R --run-id scz-all_individuals.EA-crossregion-YYYYMMDD
##
## Applies config/scz_ea_targeted.yml:reading_rule, which was fixed before any EA
## result existed, and seals the run. Every row carries gating = FALSE and
## concordance_only = TRUE: the reported quantities are an EA effect's sign
## relative to the AA effect and its nominal p. EA and AA magnitudes are never
## contrasted (AGENTS.md 7.7; config/analysis_thresholds.yml:donor_group).

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "09_schizophrenia_gwas_loci"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
man <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
smoke <- identical(toupper(man$value[man$field == "smoke_run"][1]), "TRUE")
cfg <- yaml::read_yaml(file.path(run_dir, "code", "config", "scz_ea_targeted.yml"))
pol <- donor_group_inference_policy()
regions <- names(cfg$accepted_runs)
alpha <- cfg$reading_rule$alpha

pairs <- rbindlist(lapply(regions, function(re) {
    f <- file.path(run_dir, "results", re, "pair-tests.tsv")
    if (!file.exists(f)) stop("Missing ", f, "; stage 20 did not complete for ", re)
    fread(f)
}), fill = TRUE)
targets <- fread(file.path(run_dir, "results", "targets.tsv"))
if (nrow(pairs) != nrow(targets)) stop("pair-tests rows ", nrow(pairs), " != targets ", nrow(targets))
pairs[, sign_agrees := as.logical(sign_agrees)]
pairs[, reproduces := ea_status == "tested" & sign_agrees & ea_pval < alpha]

lead <- pairs[is_lead_pair == TRUE]
n_testable_lead <- lead[ea_status == "tested", .N]
lead[, bonferroni_alpha := alpha / n_testable_lead]
lead[, reproduces_bonferroni := reproduces & ea_pval < bonferroni_alpha]
lead[, reading := fifelse(ea_status != "tested", paste0("not_testable:", ea_status),
                   fifelse(reproduces, "reproduces", fifelse(sign_agrees, "same_sign_not_nominal",
                                                             "opposite_sign")))]

by_locus <- pairs[, .(n_nominated_pairs = .N,
                      n_tested = sum(ea_status == "tested"),
                      frac_same_sign = mean(sign_agrees[ea_status == "tested"]),
                      frac_same_sign_nominal = mean(reproduces[ea_status == "tested"]),
                      ea_min_pval_same_sign = suppressWarnings(min(ea_pval[ea_status == "tested" & sign_agrees]))),
                  by = .(region, locus_id)]
loci <- fread(file.path(run_dir, "results", "loci.tsv"))
locus_table <- merge(loci, by_locus, by = c("region", "locus_id"))
locus_table <- merge(locus_table,
                     lead[, .(region, locus_id, lead_ea_status = ea_status, lead_ea_n = ea_n,
                              lead_ea_slope = ea_slope, lead_ea_slope_se = ea_slope_se,
                              lead_ea_pval = ea_pval, lead_ea_maf = ea_maf,
                              lead_sign_agrees = sign_agrees, reading, reproduces_bonferroni,
                              bonferroni_alpha)],
                     by = c("region", "locus_id"))

coloc_f <- file.path(run_dir, "results", "ea-coloc.tsv")
coloc <- if (file.exists(coloc_f)) fread(coloc_f) else data.table(status = "NOT_RUN")

label <- function(dt) {
    dt[, `:=`(gating = FALSE, concordance_only = TRUE,
              donor_group_inference = pol$donor_group_inference,
              ancestry_effect_claim_allowed = pol$ancestry_effect_claim_allowed,
              estimation_group = "EA", run_id = opts$run_id)]
    dt[]
}
pairs <- label(pairs); locus_table <- label(locus_table)

res_dir <- file.path(run_dir, "results")
write_atomic(pairs, file.path(res_dir, "ea-pair-tests.tsv"))
write_atomic(locus_table, file.path(res_dir, "ea-locus-reading.tsv"))
summary_dt <- data.table(
    n_loci = nrow(locus_table),
    n_lead_testable = n_testable_lead,
    n_lead_reproduce_nominal = lead[reading == "reproduces", .N],
    n_lead_reproduce_bonferroni = lead[reproduces_bonferroni == TRUE, .N],
    n_lead_same_sign = lead[ea_status == "tested" & sign_agrees == TRUE, .N],
    n_pairs = nrow(pairs), n_pairs_tested = pairs[ea_status == "tested", .N],
    n_pairs_same_sign = pairs[ea_status == "tested" & sign_agrees == TRUE, .N],
    n_coloc_run = if ("status" %in% names(coloc)) coloc[status == "OK", .N] else 0L,
    n_coloc_pp4_threshold = if ("pp4_meets_threshold" %in% names(coloc)) coloc[pp4_meets_threshold == TRUE, .N] else 0L)
write_atomic(summary_dt, file.path(res_dir, "ea-summary.tsv"))
writeLines(c(
    "Interpretation constraints carried by this run:",
    "  - NOT A GATE. Module 09's locked decisions are unchanged; coloc here is",
    "    exploratory support and never part of scz_application_retention.",
    "  - CONCORDANCE ONLY (AGENTS.md 7.7). An EA effect is read only for its sign",
    "    relative to the AA effect and its nominal p. EA and AA magnitudes are not",
    "    contrasted, and a non-replication is not evidence of an ancestry-specific",
    "    effect: EA cells are 129 / 55 / 60 donors.",
    "  - The model is the locked M3a, with snpPC1-5 and methPC1-5 re-estimated",
    "    inside the EA donors.",
    "  - Loci are the prespecified illustrative loci of the accepted AA runs; pairs",
    "    are those significant in the AA primary. Nothing was selected on EA data,",
    "    except coloc eligibility, by a rule fixed in config.",
    "  - Coloc is PGC3 European x EA meQTL, the LD-matched pairing. coloc.abf assumes",
    "    one causal variant per trait; small EA n makes a low PP4 weak evidence.",
    "  - GWAS loci are defined from European-ancestry summary statistics (AGENTS.md",
    "    7.8 rule 2)."
), file.path(res_dir, "interpretation-constraints.txt"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "session-info.txt"))
print(locus_table[, .(region, rank, locus_id, index_snp_prioritized, lead_cpg_id, lead_aa_pval,
                      lead_ea_n, lead_ea_pval, lead_sign_agrees, reading, n_tested,
                      frac_same_sign, frac_same_sign_nominal)])
print(summary_dt)
print(coloc)

if (!smoke) {
    ## The per-region working files (dosages, PCA, staged genotypes) are inputs
    ## to the results, not results; they are removed before sealing so the run
    ## carries no individual-level genotype dosages.
    unlink(file.path(run_dir, "work"), recursive = TRUE)
    append_manifest(list(dir = run_dir), as.list(c(
        lapply(summary_dt, as.character),
        git_commit = git_commit(), git_dirty = as.character(git_dirty()),
        sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))))
    close_run(list(dir = run_dir))
    message("[09.22] sealed ", opts$run_id)
} else {
    message("[09.22] smoke run, not sealed")
}
