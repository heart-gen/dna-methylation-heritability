#!/usr/bin/env Rscript
#### 02c -- reconcile, then the existence and ordering readouts ####
##
## Every (VMR, mode) unit must be completed, QC-failed with a recorded reason,
## or failed; reconcile() stops on anything unexplained (AGENTS.md 9).
##
## Then, per REML mode, on Module 02-eligible VMRs whose fit converged in that
## mode -- each mode on its own converged set, never one mode filling another's
## gaps:
##   existence  mean cis-GREML estimate across VMRs, delete-one-chromosome
##              jackknife CI. 02b found the unconstrained estimator unbiased on
##              average, so a mean above zero is evidence that local genetic
##              control exists in these data; no single VMR is called.
##   ordering   Spearman(cis-GREML estimate, Module 02 score), jackknife CI; the
##              same against Module 02's HE and BSLMM features, descriptively.
##   profile    mean estimate per Module 02 score decile, jackknife SE.
##   selection  convergence rate per score decile, so a comparison restricted
##              to converged fits shows what it conditioned on.
##   LRT        fraction with p < 0.05, DESCRIPTIVE ONLY: 02b measured CI
##              coverage ~0.75 at true h2 = 0, so it is not a calibrated test
##              per VMR and nothing is classified on it.
##
## Usage: Rscript _h/02_summarize.R --run-id <id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02c_cis_greml_sensitivity", "_h")),
                 "cgs_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
run <- list(run_id = opts$run_id, dir = run_dir)
cfg <- load_cgs_config(run_dir)
man <- read_run_manifest(run_dir)
smoke <- identical(man[["smoke_run"]], "TRUE")
tasks <- fread(file.path(run_dir, "config", "vmr-tasks.tsv"))
modes <- fread(file.path(run_dir, "config", "reml-modes.tsv"))
expected <- CJ(vmr_id = tasks$vmr_id, reml_mode = modes$reml_mode)
expected[, unit := unit_id(vmr_id, reml_mode)]
st <- rbindlist(lapply(list.files(file.path(run_dir, "status"), full.names = TRUE),
                       fread, colClasses = "character"), fill = TRUE)
reconcile(expected$unit, completed = st[status == "completed", unit],
          qc_failed = st[status == "qc_failed", unit],
          failed = st[status == "failed", unit], run = run, allow_failures = smoke)

fits <- rbindlist(lapply(list.files(file.path(run_dir, "results"), "^fits_task_",
                                    full.names = TRUE), fread), fill = TRUE)
loci <- rbindlist(lapply(list.files(file.path(run_dir, "results"), "^loci_task_",
                                    full.names = TRUE), fread), fill = TRUE)
m02 <- tasks[, .(vmr_id, module02_eligible, module02_exclusion_reason,
                 local_snp_contribution_score, local_snp_contribution_score_z,
                 he_h2, bslmm_pve)]
vmr <- merge(fits, m02, by = "vmr_id", all.x = TRUE)
vmr[, chrom_block := sub("^chr", "", chrom)]
score_col <- cfg$summaries$score_column
ndec <- as.integer(cfg$summaries$score_deciles)
vmr[, score_decile := cut(get(score_col), breaks = seq(0, 1, length.out = ndec + 1),
                          labels = FALSE, include.lowest = TRUE)]
cgs_flags(vmr)
write_atomic(vmr, file.path(run_dir, "summary", "vmr-cis-greml.tsv"))
write_atomic(loci, file.path(run_dir, "summary", "locus-snp-check.tsv"))

spear <- function(a, b) function(d) suppressWarnings(
    stats::cor(d[[a]], d[[b]], method = "spearman", use = "complete.obs"))
out_exist <- list(); out_order <- list(); out_prof <- list(); out_sel <- list()
for (m in modes$reml_mode) {
    el <- vmr[reml_mode == m & module02_eligible %in% c(TRUE, "TRUE")]
    cv <- el[converged == TRUE & is.finite(h2_hat)]
    role <- modes[reml_mode == m, role]
    ex <- block_jackknife_stat(cv, "chrom_block", function(d) mean(d$h2_hat))
    out_exist[[m]] <- cbind(data.table(reml_mode = m, reml_role = role,
        statistic = "mean_cis_greml_h2", n_eligible = nrow(el), n_converged = nrow(cv),
        convergence_rate = nrow(cv) / max(1, nrow(el)),
        median_h2_hat = stats::median(cv$h2_hat),
        frac_h2_hat_positive = mean(cv$h2_hat > 0),
        frac_lrt_p_below_0_05_descriptive = mean(cv$pval < 0.05, na.rm = TRUE),
        lrt_used_for_classification = FALSE), ex)
    cols <- c(score_col, unlist(cfg$summaries$module02_estimator_columns))
    out_order[[m]] <- rbindlist(lapply(cols, function(cc) {
        d <- cv[is.finite(get(cc))]
        cbind(data.table(reml_mode = m, reml_role = role, against = cc,
                         role_of_comparison = if (cc == score_col) "primary_ordering"
                                              else "descriptive_estimator_concordance",
                         n = nrow(d)),
              block_jackknife_stat(d, "chrom_block", spear("h2_hat", cc)))
    }))
    out_prof[[m]] <- cv[, cbind(data.table(n = .N),
                                block_jackknife_stat(.SD, "chrom_block",
                                                     function(d) mean(d$h2_hat))),
                        by = score_decile][order(score_decile)][, `:=`(reml_mode = m,
                                                                       reml_role = role)]
    out_sel[[m]] <- el[, .(n_eligible = .N, n_converged = sum(converged == TRUE),
                           convergence_rate = mean(converged == TRUE)),
                       by = score_decile][order(score_decile)][, reml_mode := m]
}
ex <- rbindlist(out_exist); od <- rbindlist(out_order)
pf <- rbindlist(out_prof); sl <- rbindlist(out_sel)
for (d in list(ex, od, pf, sl)) { cgs_flags(d); d[, region := man[["region"]]] }
write_atomic(ex, file.path(run_dir, "summary", "existence.tsv"))
write_atomic(od, file.path(run_dir, "summary", "ordering.tsv"))
write_atomic(pf, file.path(run_dir, "summary", "decile-profile.tsv"))
write_atomic(sl, file.path(run_dir, "summary", "convergence-by-decile.tsv"))
append_manifest(run, list(summarized_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
print(ex[, .(reml_mode, n_eligible, n_converged, estimate, ci_low, ci_high)])
print(od[, .(reml_mode, against, n, estimate, ci_low, ci_high)])
