#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- locus-level transcriptional coupling ####
##
## Usage:
##   Rscript _h/06_transcriptional_coupling.R --run-id scz-AA-caudate-YYYYMMDD
##
## AGENTS.md 7.8 retention criterion: "at least one locus with transcriptional
## coupling". Module 07 already tested coupling for every VMR of this cell and
## its run is accepted; this stage projects that result onto loci. It runs no
## new association, which is why it cannot inflate the coupling evidence.
##
## Module 07's accepted result is heterogeneous by design: nearest-gene
## expression coupling in all three regions, PSI strong in caudate, thin in
## DLPFC (24 VMRs) and null in hippocampus, ABC underpowered, and the internal
## LIBD eQTL arm switched off pending its QC repair. Per-modality counts are
## therefore carried through to the locus level rather than collapsed, so a
## locus supported only by a thin modality is visible as such.
##
## Coupling is an association between genetic regulation and a local
## transcriptional signal. It is NOT mediation and implies no causal ordering;
## that constraint travels with the table.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09_schizophrenia_risk_application", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

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

links <- fread(file.path(run_dir, "results", "locus-vmr-links.tsv"))
tsc_run <- mf("upstream_tx_coupling_run_id")
tsc_dir <- file.path(repo_root(), "07_transcription_splicing_coupling", "_m",
                     "runs", tsc_run, "results")
summ_files <- list.files(tsc_dir, pattern = "-vmr-summary\\.tsv$", full.names = TRUE)
if (length(summ_files) == 0) stop("No Module 07 vmr-summary tables in ", tsc_dir)
coupling <- rbindlist(lapply(summ_files, fread), use.names = TRUE, fill = TRUE)
coupling[, any_sig_fdr := as.logical(any_sig_fdr)]

## Which VMRs of this locus were tested by Module 07, and which were coupled.
vmr_locus <- unique(links[, .(locus_id, vmr_id)])
joined <- merge(vmr_locus, coupling, by = "vmr_id", allow.cartesian = TRUE)

if (nrow(joined) == 0) {
    write_atomic(data.table(), file.path(run_dir, "results", "locus-coupling.tsv"))
    message("[09] no linked VMR was tested by Module 07; coupling table empty")
    quit(status = 0)
}

per_modality <- joined[, .(
    n_vmrs_tested = uniqueN(vmr_id),
    n_vmrs_coupled = uniqueN(vmr_id[any_sig_fdr]),
    min_fdr = min(min_fdr, na.rm = TRUE),
    max_abs_r = max(max_abs_r, na.rm = TRUE)
), by = .(locus_id, modality)]

per_locus <- joined[, .(
    n_vmrs_tested_any_modality = uniqueN(vmr_id),
    n_vmrs_coupled_any_modality = uniqueN(vmr_id[any_sig_fdr]),
    modalities_tested = paste(sort(unique(modality)), collapse = ","),
    modalities_coupled = paste(sort(unique(modality[any_sig_fdr])), collapse = ","),
    min_fdr_any_modality = min(min_fdr, na.rm = TRUE),
    has_transcriptional_coupling = any(any_sig_fdr)
), by = locus_id]

out <- merge(per_locus, dcast(per_modality, locus_id ~ modality,
                              value.var = "n_vmrs_coupled", fill = 0),
             by = "locus_id", all.x = TRUE)
out[, `:=`(run_id = opts$run_id, cohort = cohort, region = region,
           upstream_tx_coupling_run_id = tsc_run,
           permitted_claim =
               "the locus contains a VMR with a local transcriptional association",
           forbidden_claim =
               "methylation mediates the genetic effect on expression or splicing")]
setorder(out, -has_transcriptional_coupling, min_fdr_any_modality, locus_id)
write_atomic(out, file.path(run_dir, "results", "locus-coupling.tsv"))
write_atomic(per_modality,
             file.path(run_dir, "results", "locus-coupling-by-modality.tsv"))

print(out[, .(n_loci = .N,
              n_loci_coupled = sum(has_transcriptional_coupling))])
message("[09] ", sum(out$has_transcriptional_coupling), "/", nrow(out),
        " loci carry transcriptional coupling")
