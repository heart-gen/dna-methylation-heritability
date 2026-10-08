#!/usr/bin/env Rscript
#### 06_partitioned_heritability -- external-annotation benchmark: combine and seal ####
##
## Usage:
##   Rscript _h/13_external_summarize.R --run-id sldsc-AA-external-YYYYMMDD
##
## Combines stage 12's per-(annotation, trait) metrics and refuses a partial
## grid: every staged annotation x every frozen trait must have completed with a
## finite, positive tau SE. Beside them it places the accepted two-annotation
## VMR_TESTED rows, unchanged, so the like-for-like standalone refit and the
## reported model sit in one table.
##
## Nothing here is a gate on any accepted result. q-values are BH across the 8
## traits WITHIN each annotation and are descriptive.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "06_partitioned_heritability"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
man <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
smoke <- identical(toupper(man$value[man$field == "smoke_run"][1]), "TRUE")

ph <- yaml::read_yaml(file.path(run_dir, "code", "config", "partitioned_heritability.yml"))
ext <- yaml::read_yaml(file.path(run_dir, "code", "config", "sldsc_external_annotations.yml"))
traits <- data.table(trait = vapply(ph$traits, `[[`, character(1), "name"),
                     trait_class = vapply(ph$traits, `[[`, character(1), "class"),
                     trait_label = vapply(ph$traits, `[[`, character(1), "label"))
listing <- fread(file.path(run_dir, "annotation", "annotation-list.tsv"))
ready <- listing[status == "ready"]

## A smoke run may restrict the trait fan-out (SMOKE_TRAITS); the grid it is
## held to is then the traits it submitted. Production is held to all of them.
if (smoke) {
    jobs <- fread(file.path(run_dir, "submitted-jobs.tsv"))
    submitted <- unique(sub(".*:", "", jobs[step == 12, script]))
    traits <- traits[trait %in% submitted]
}
grid <- CJ(annotation = ready$annotation, trait = traits$trait)
metrics <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
    f <- file.path(run_dir, "results", "sldsc", grid$annotation[i],
                   paste0(grid$trait[i], ".metrics.tsv"))
    if (!file.exists(f)) return(NULL)
    fread(f)
}), fill = TRUE)
missing <- grid[!metrics, on = c("annotation", "trait")]
if (nrow(missing)) {
    print(missing)
    stop(nrow(missing), " annotation x trait regressions did not complete; a partial ",
         "grid is refused.")
}
bad_se <- metrics[!(is.finite(tau_se) & tau_se > 0)]
if (nrow(bad_se)) {
    print(bad_se[, .(annotation, trait, tau, tau_se)])
    stop("Non-finite or non-positive tau SE; the regression did not identify the model.")
}
metrics <- merge(metrics, traits, by = "trait")
metrics <- merge(metrics, listing[, .(annotation, role, label)], by = "annotation")
metrics[, q_bh_within_annotation := p.adjust(enrichment_p, method = ext$descriptive_fdr_method),
        by = annotation]
metrics[, `:=`(model = "baselineLD_v2.2_plus_one_binary_annotation",
               ld_reference_arm = ph$ld_reference_arm, gating = FALSE,
               total_h2_z = total_h2 / total_h2_se, run_id = opts$run_id)]

## The accepted two-annotation rows, read unchanged.
accepted <- rbindlist(lapply(names(ext$accepted_vmr_runs), function(re) {
    f <- file.path(repo_root(), MODULE, "_m", "runs", ext$accepted_vmr_runs[[re]],
                   "results", "sldsc-membership-metrics.tsv")
    d <- fread(f)
    d[, .(annotation = paste0("VMR_TESTED_", toupper(re), "__ACCEPTED_TWO_ANNOTATION_MODEL"),
          trait, prop_snps, prop_h2, prop_h2_se, enrichment, enrichment_se, enrichment_p,
          tau, tau_se, tau_z, tau_p_two_sided = tau_p_nominal, total_h2, total_h2_se,
          role = "vmr_membership_accepted_joint_model",
          model = "baselineLD_v2.2_plus_VMR_TESTED_plus_LOCAL_SNP_CONTRIBUTION_Z",
          source_run_id = ext$accepted_vmr_runs[[re]])]
}))

res_dir <- file.path(run_dir, "results")
write_atomic(metrics, file.path(res_dir, "external-annotation-metrics.tsv"))
write_atomic(accepted, file.path(res_dir, "vmr-membership-accepted-metrics.tsv"))
write_atomic(listing, file.path(res_dir, "annotation-coverage.tsv"))

compact <- metrics[, .(annotation, role, trait, trait_class,
                       prop_snps = signif(prop_snps, 3), prop_h2 = signif(prop_h2, 3),
                       enrichment = signif(enrichment, 3), enrichment_se = signif(enrichment_se, 3),
                       enrichment_p = signif(enrichment_p, 3), tau_star = signif(tau_star, 3),
                       tau_star_se = signif(tau_star_se, 3), tau_z = signif(tau_z, 3),
                       q_bh_within_annotation = signif(q_bh_within_annotation, 3))]
setorder(compact, annotation, trait)
write_atomic(compact, file.path(res_dir, "external-annotation-summary.tsv"))

skipped <- listing[status != "ready"]
writeLines(c(
    "Interpretation constraints carried by this run:",
    "  - NON-GATING. The accepted Module 06 decision (sldsc_supports_brain_enrichment",
    "    = FALSE, sldsc-AA-*-20260925) is unchanged. This run is the positive control",
    "    that run lacked: same baseline model, LD reference, weights, frq and munged",
    "    GWAS files; each annotation enters alone on top of baselineLD v2.2.",
    "  - Reading: if a Rizzardi NeuN+ CG-DMR annotation of comparable footprint is",
    "    enriched for the brain traits and the standalone VMR_TESTED is not, the",
    "    pipeline can detect enrichment at this footprint and the VMR null is not a",
    "    pure power null. If neither is, the VMR null stays 'no detectable",
    "    enrichment at this footprint'.",
    "  - Enrichment is interpretable for every annotation here: all are binary.",
    "  - q-values are BH across the 8 traits within an annotation, descriptive only.",
    "  - EUR LD and EUR-dominated GWAS; the annotations are genomic features, not",
    "    donor properties (Module 06 README, LD reference and ancestry).",
    if (nrow(skipped)) paste0("  - NOT COVERED: ", paste0(skipped$annotation, " (", skipped$status, ")",
                                                         collapse = "; ")) else NULL
), file.path(res_dir, "interpretation-constraints.txt"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "session-info.txt"))

print(compact[trait_class == "brain" | trait %in% c("asthma", "cad")], nrows = 200)
if (!smoke) {
    append_manifest(list(dir = run_dir), list(
        n_annotations_ready = nrow(ready),
        n_annotations_not_covered = nrow(skipped),
        annotations_not_covered = paste(skipped$annotation, collapse = ","),
        n_regressions = nrow(metrics),
        git_commit = git_commit(),
        git_dirty = as.character(git_dirty()),
        sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
    close_run(list(dir = run_dir))
    message("[06x] sealed ", opts$run_id)
} else {
    message("[06x] smoke run, not sealed")
}
