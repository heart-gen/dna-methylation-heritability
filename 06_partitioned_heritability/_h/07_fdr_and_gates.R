#!/usr/bin/env Rscript
#### 06_partitioned_heritability -- pool traits, apply FDR, decide ####
##
## Usage:
##   Rscript _h/07_fdr_and_gates.R --run-id sldsc-AA-caudate-YYYYMMDD
##
## The legacy fdr_correction.py hardcoded `allowed_diseases = {"smoking","ad",
## "pd","scz"}` and globbed whatever .results files happened to exist, so the
## FDR family was whichever subset of runs had completed. Here the family is the
## frozen trait list from config: every declared trait must have produced
## metrics, or the run fails rather than correcting across a smaller family and
## reporting optimistic q-values.
##
## Stage 06 now emits two rows per trait, one per annotation in the
## two-annotation S-LDSC model (see _h/annotations.py). PRESPECIFIED, before any
## two-annotation result was seen: the FDR family stays exactly the frozen
## traits on the SCORE annotation, which carries the module's hypothesis. The
## VMR_TESTED membership tau is a different question -- are tested VMRs enriched
## at all -- and is written to its own descriptive table with no q-value. That
## keeps the family size at the declared trait count and leaves every q-value
## the primary hypothesis depends on computed over the same family as before.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE", file.path(Sys.getenv("V2_REPO_ROOT", "."), "06_partitioned_heritability", "_h")), "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "06_partitioned_heritability"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
## Select outside the data.table `[` -- see the note in 01_build_annotation.R.
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
smoke <- identical(mf("smoke_run"), "TRUE")

ph <- load_run_config("partitioned_heritability", run_dir)
declared <- vapply(ph$traits, `[[`, character(1), "name")
trait_class <- setNames(vapply(ph$traits, `[[`, character(1), "class"), declared)
trait_label <- setNames(vapply(ph$traits, `[[`, character(1), "label"), declared)

metric_dir <- file.path(run_dir, "results", "sldsc")
files <- file.path(metric_dir, paste0(declared, ".metrics.tsv"))
present <- file.exists(files)

if (!all(present) && isTRUE(ph$gates$require_all_traits_completed) && !smoke) {
    stop("Missing S-LDSC metrics for: ",
         paste(declared[!present], collapse = ", "),
         "\n  The frozen trait list defines the FDR family; correcting across ",
         "a partial family would understate every q-value.")
}

all_rows <- rbindlist(lapply(files[present], fread), fill = TRUE)
if (nrow(all_rows) == 0) stop("No S-LDSC metrics found in ", metric_dir)

## Key on the ROLE, not the annotation name: stage 06 stamps the role on every
## row, so the annotation's name lives in exactly one place (_h/annotations.py)
## and cannot drift between the Python and R halves of the module.
if (!all(c("annotation_role", "annotation_name") %in% names(all_rows))) {
    stop("Stage 06 metrics carry no `annotation_role`/`annotation_name`. ",
         "These are ",
         "one-annotation metrics from the superseded single-column model, ",
         "whose tau is not conditional on VMR membership and is therefore not ",
         "interpretable as the within-VMR gradient this module reports. ",
         "Rerun stages 03, 05 and 06.")
}

REQUIRED_ROLES <- c("membership", "score")
role_counts <- all_rows[, .N, by = .(trait, annotation_role)]
missing_roles <- setdiff(REQUIRED_ROLES, unique(all_rows$annotation_role))
dup_roles <- role_counts[N > 1]

all_rows[, trait_class := trait_class[trait]]
all_rows[, trait_label := trait_label[trait]]

## The FDR family. The membership rows are split off BEFORE p.adjust() so the
## family cannot silently double in size.
res <- all_rows[annotation_role == "score"]
member <- all_rows[annotation_role == "membership"]
if (nrow(res) == 0) {
    stop("No rows with annotation_role == 'score' in ", metric_dir,
         "; the primary annotation produced no result.")
}

## Two-tailed p from the tau z-score. The coefficient z is the statistic that is
## conditional on baselineLD, so it -- not the enrichment ratio -- is what the
## FDR family is built on. Enrichment and its own p are reported alongside, as
## the legacy interpreting_sldsc_results.md requires.
res[, tau_p := 2 * pnorm(-abs(as.numeric(tau_z)))]

## config names the procedure (`fdr_bh`), p.adjust names an implementation
## (`BH`). Translate here rather than writing p.adjust's spelling into the
## locked config: the config should say which multiple-testing procedure was
## prespecified, not which R function argument realizes it. Unrecognized values
## abort -- defaulting to BH would let a typo in a pi_locked config silently
## change the correction that every reported q-value depends on.
FDR_METHOD <- c(fdr_bh = "BH", fdr_by = "BY", bonferroni = "bonferroni",
                holm = "holm", none = "none")
if (!isTRUE(ph$fdr_method %in% names(FDR_METHOD))) {
    stop("config fdr_method is ", deparse(ph$fdr_method),
         "; expected one of: ", paste(names(FDR_METHOD), collapse = ", "))
}
padjust_method <- unname(FDR_METHOD[[ph$fdr_method]])

res[, tau_q := p.adjust(tau_p, method = padjust_method)]
res[, enrichment_q := p.adjust(as.numeric(enrichment_p), method = padjust_method)]
res[, in_fdr_family := TRUE]
res[, fdr_family := ph$fdr_family]
res[, fdr_method := ph$fdr_method]
res[, n_tests_in_family := .N]

stamp <- function(dt) {
    dt[, cohort := mf("cohort")]
    dt[, region := mf("region")]
    dt[, ld_reference_arm := mf("ld_reference_arm")]
    dt[, run_id := opts$run_id]
    dt[]
}
stamp(res)

setorder(res, tau_p)
write_atomic(res, file.path(run_dir, "results", "sldsc-metrics.tsv"))

## The membership annotation, reported and NOT corrected. Its nominal p is
## computed because a reader will want it, but it carries no q-value: adding it
## to the frozen family would change the family size, and reporting a q from a
## family of one would be worse than reporting none. Kept in a separate file so
## no consumer of sldsc-metrics.tsv -- Module 09 included -- can pick these rows
## up as if they were the primary hypothesis.
if (nrow(member) > 0) {
    member[, tau_p_nominal := 2 * pnorm(-abs(as.numeric(tau_z)))]
    member[, in_fdr_family := FALSE]
    member[, fdr_family := NA_character_]
    member[, interpretation := paste(
        "descriptive: are tested VMRs enriched at all.",
        "Not in the frozen FDR family; nominal p only.")]
    stamp(member)
    setorder(member, tau_p_nominal)
    write_atomic(member, file.path(run_dir, "results",
                                   "sldsc-membership-metrics.tsv"))
}

## ------------------------------------------------------------------- gates
fail <- character()

if (!all(present) && !smoke) fail <- c(fail, "INCOMPLETE_TRAIT_FAMILY")

## The estimand gate. A run whose model lacks the membership annotation produces
## a tau that conflates "not a tested VMR" with "a tested VMR of median local
## genetic control", so it cannot be read as the within-VMR gradient this module
## reports -- null or otherwise. That is a QC failure and not a caveat.
if (length(missing_roles) > 0) {
    fail <- c(fail, paste0("ANNOTATION_ROLE_MISSING(",
                           paste(missing_roles, collapse = ","), ")"))
}
if (nrow(dup_roles) > 0) {
    fail <- c(fail, paste0("DUPLICATE_ANNOTATION_ROW(",
                           paste(paste0(dup_roles$trait, ":",
                                        dup_roles$annotation_role),
                                 collapse = ","), ")"))
}
incomplete_annot <- setdiff(res$trait, member$trait)
if (length(incomplete_annot) > 0) {
    fail <- c(fail, paste0("TRAIT_MISSING_MEMBERSHIP_ROW(",
                           paste(incomplete_annot, collapse = ","), ")"))
}

if (isTRUE(ph$gates$require_finite_tau_se)) {
    ## Both annotations: they are fitted jointly, so a degenerate SE on either
    ## means the joint fit did not identify the model, and the score's tau is
    ## then no longer the conditional quantity it is reported as.
    bad <- all_rows[!is.finite(as.numeric(tau_se)) | as.numeric(tau_se) <= 0]
    if (nrow(bad) > 0) {
        fail <- c(fail, paste0("NONFINITE_TAU_SE(",
                               paste(paste0(bad$trait, ":",
                                            bad$annotation_role),
                                     collapse = ","), ")"))
    }
}

## A trait whose total h2 is not distinguishable from zero carries no
## heritability to partition, so its enrichment is uninterpretable rather than
## null. Such traits are dropped from interpretation, not from the family.
res[, total_h2_z := as.numeric(total_h2) / as.numeric(total_h2_se)]
res[, h2_interpretable := is.finite(total_h2_z) &
        total_h2_z >= ph$gates$min_total_h2_z]
if (!any(res$h2_interpretable)) {
    fail <- c(fail, "NO_TRAIT_WITH_INTERPRETABLE_TOTAL_H2")
}

ann <- fread(file.path(run_dir, "annotation", "annotation-summary.tsv"))
if (!isTRUE(ann$annotation_continuous[1]) ||
    isTRUE(ann$annotation_thresholded[1]) ||
    isTRUE(ann$absolute_pve_used[1])) {
    fail <- c(fail, "ANNOTATION_NOT_CONTINUOUS_OR_USES_BANNED_QUANTITY")
}

decision <- if (length(fail) > 0) {
    paste0("FAIL_PARTITIONED_H2_QC:", paste(fail, collapse = ";"))
} else if (smoke) {
    "PASS_SMOKE_ONLY_NOT_ACCEPTABLE"
} else {
    "PASS_PARTITIONED_H2_QC"
}

brain_sig <- res[trait_class == "brain" & h2_interpretable &
                 tau_q < ph$fdr_threshold]
control_sig <- res[trait_class == "control" & h2_interpretable &
                   tau_q < ph$fdr_threshold]

dec <- data.table(
    run_id = opts$run_id, cohort = mf("cohort"), region = mf("region"),
    decision = decision,
    smoke_run = smoke,
    ld_reference_arm = mf("ld_reference_arm"),
    n_traits_declared = length(declared),
    n_traits_completed = sum(present),
    n_traits_h2_interpretable = sum(res$h2_interpretable),
    n_brain_significant = nrow(brain_sig),
    n_control_significant = nrow(control_sig),
    brain_significant = paste(brain_sig$trait, collapse = ","),
    control_significant = paste(control_sig$trait, collapse = ","),
    fdr_family = ph$fdr_family, fdr_threshold = ph$fdr_threshold,
    ## The estimand, recorded on the decision so a reader of this row alone can
    ## tell whether tau was conditional on VMR membership.
    annotations_in_model = paste(sort(unique(all_rows$annotation_name)),
                                 collapse = ","),
    fdr_family_annotation = paste(sort(unique(res$annotation_name)),
                                  collapse = ","),
    tau_conditional_on_vmr_membership = "membership" %in%
        unique(all_rows$annotation_role),
    failures = paste(fail, collapse = ";"),
    absolute_pve_interpretation_allowed = FALSE,
    ## AGENTS.md 7.8 / config/analysis_thresholds.yml: Module 09 needs to know
    ## whether S-LDSC added anything, and that is a property of THIS run.
    sldsc_supports_brain_enrichment = nrow(brain_sig) > 0 && !smoke &&
        length(fail) == 0
)
write_atomic(dec, file.path(run_dir, "results", "partitioned-h2-decision.tsv"))

message("[06] primary annotation (in the frozen FDR family):")
print(res[, .(trait, trait_class, enrichment, enrichment_p, tau_z, tau_q,
              h2_interpretable)])
if (nrow(member) > 0) {
    message("[06] membership annotation (descriptive, not FDR-corrected):")
    print(member[, .(trait, trait_class, m_5_50, enrichment, tau_z,
                     tau_p_nominal)])
}
message("[06] decision ", decision)
if (startsWith(decision, "FAIL")) {
    quit(status = 1)
}
