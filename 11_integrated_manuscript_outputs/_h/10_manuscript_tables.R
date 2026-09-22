#### 11 / Manuscript tables, registry and provenance products ####
##
## AGENTS.md 7.11 makes this module the one source of truth and names its
## products. Table 1 and the cohort QC come from 04/05; everything else in that
## list is built here:
##
##   tables/manuscript-number-registry.tsv   every citable number -> its panel,
##                                           run, table, column and filter
##   tables/analysis-to-claim-matrix.tsv     claim -> module -> accepted run ->
##                                           decision token -> permitted claim
##   tables/exclusions-and-denominators.tsv  donors, VMRs called, eligible, why
##   tables/supplementary-table-index.tsv    the _m/combined/ deliverables
##   tables/software-and-run-manifest.tsv    environment and upstream run IDs
##   tables/analysis-to-claim-matrix.tex     booktabs fragment
##
## MUST RUN AFTER THE FIGURES. The registry is built by reading the
## source_data/ tables the figure builders wrote, which is what makes it a
## record of what was actually rendered rather than a second, driftable list.
## It is also Supplementary Data 14.
##
## Usage:
##   Rscript 10_manuscript_tables.R --cohort AA --run-id fig-all-YYYYMMDD

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "11_integrated_manuscript_outputs", "_h",
                 "00_figure_theme.R"))

opts <- parse_v2_args(require = c("cohort", "run_id"))
cohort  <- opts$cohort
cohorts <- load_config("cohorts")
regions <- cohorts$regions
arms    <- cohorts$arms

module_root <- file.path(V2_ROOT, "11_integrated_manuscript_outputs")
run_dir  <- if (!is.null(opts$out_dir)) opts$out_dir else
    file.path(module_root, "_m", "runs", opts$run_id)
fig_dir   <- file.path(run_dir, "figures")
data_dir  <- file.path(run_dir, "source_data")
table_dir <- file.path(run_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

SCRIPT <- "11_integrated_manuscript_outputs/_h/10_manuscript_tables.R"

## ================================================ 1. number registry (SD14)
srcs <- list.files(data_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(srcs) == 0) {
    stop("No source_data/ tables in ", run_dir, ". The figure builders must ",
         "run before 10_manuscript_tables.R -- the registry is built from what ",
         "they actually rendered.")
}

registry <- rbindlist(lapply(srcs, function(f) {
    d <- fread(f)
    stem <- sub("\\.tsv$", "", basename(f))
    ## "<figure stem>_panel_<tag>" -> figure and panel. Anything without the
    ## panel marker is a whole-figure companion table (a decision row, a claim
    ## table) and is recorded with panel = NA rather than guessed at.
    panel <- sub("^.*_panel_", "", stem)
    has_panel <- grepl("_panel_", stem)
    data.table(
        figure        = if (has_panel) sub("_panel_.*$", "", stem) else stem,
        panel         = if (has_panel) panel else NA_character_,
        source_file   = paste0("source_data/", basename(f)),
        n_rows        = nrow(d),
        n_columns     = ncol(d),
        source_run_id = if ("source_run_id" %in% names(d)) d$source_run_id[1] else NA_character_,
        source_table  = if ("source_table"  %in% names(d)) d$source_table[1]  else NA_character_,
        source_script = if ("source_script" %in% names(d)) d$source_script[1] else NA_character_,
        row_filter    = if ("row_filter"    %in% names(d)) d$row_filter[1]    else NA_character_,
        columns       = paste(setdiff(names(d), c("source_run_id", "source_table",
                                                  "source_script", "row_filter")),
                              collapse = ";"))
}))
setorder(registry, figure, panel, na.last = TRUE)

## Every rendered figure must appear. close_run() enforces this too, but the
## registry is what the manuscript cites, so it checks independently.
rendered <- unique(sub("\\.pdf$", "",
                       list.files(fig_dir, pattern = "\\.pdf$")))
missing <- setdiff(rendered, unique(registry$figure))
if (length(missing) > 0) {
    stop("Figures with no source data, so no citable numbers: ",
         paste(missing, collapse = ", "))
}
write_atomic(registry, file.path(table_dir, "manuscript-number-registry.tsv"))
message("[table] manuscript-number-registry.tsv (", nrow(registry), " panels across ",
        length(unique(registry$figure)), " figures)")

## ============================================= 2. analysis-to-claim matrix
##
## Sourced from each module's README acceptance table -- the record of
## decision that 00_shared/gates.R parses -- rather than retyped here, so the
## matrix cannot drift from the gate.
MODULES <- c("01_vmr_catalog", "01b_estimation_cells", "02_local_genetic_variance",
             "03_local_snp_prediction", "04_repeat_repressive_architecture",
             "05_cpg_meqtl_burden", "06_partitioned_heritability",
             "07_transcription_splicing_coupling", "08_region_donor_generalization",
             "09_schizophrenia_risk_application", "09b_aging_application",
             "10_environmental_exploratory")

claims <- rbindlist(lapply(MODULES, function(m) {
    a <- tryCatch(read_accepted_runs(m), error = function(e) data.table())
    if (nrow(a) == 0) {
        return(data.table(module = m, run_id = NA_character_, cohort = NA_character_,
                          region = NA_character_, decision = "NO ACCEPTED RUN",
                          accepted_on = NA_character_, notes = NA_character_))
    }
    keep <- intersect(c("run_id", "cohort", "region", "decision", "accepted_on",
                        "notes"), names(a))
    d <- a[, ..keep]
    d[, module := m][]
}), fill = TRUE)
setcolorder(claims, "module")

## Module 04 additionally emits the permitted-claim string per outcome, which
## is the sentence the manuscript may write. It lives on the gate host.
rra <- tryCatch({
    host <- require_accepted_upstream("04_repeat_repressive_architecture",
                                      cohort, "caudate")$run_id
    f <- file.path(V2_ROOT, "04_repeat_repressive_architecture", "_m", "runs",
                   host, "results", "interpretation-claims.tsv")
    d <- fread(f); d[, module := "04_repeat_repressive_architecture"][]
}, error = function(e) data.table())

write_atomic(claims, file.path(table_dir, "analysis-to-claim-matrix.tsv"))
if (nrow(rra) > 0) {
    write_atomic(rra, file.path(table_dir, "permitted-claims-repeat-chromatin.tsv"))
}
message("[table] analysis-to-claim-matrix.tsv (", nrow(claims), " accepted runs across ",
        length(MODULES), " modules)")

## ============================================ 3. exclusions and denominators
##
## AGENTS.md 11: always report denominators, exclusions, region, donor group,
## VMR set and the exact metric. One table, so a reviewer does not have to
## reassemble it from five figure legends.
denom <- rbindlist(lapply(arms, function(a) rbindlist(lapply(regions, function(r) {
    cat_run <- tryCatch(require_accepted_upstream("01_vmr_catalog", a, r),
                        error = function(e) NULL)
    if (is.null(cat_run)) return(NULL)
    vmr <- fread(file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs",
                           cat_run$run_id, "vmr", "vmr_catalog.tsv"))
    don <- fread(file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs",
                           cat_run$run_id, "vmr", "donors_plink.txt"),
                 header = FALSE)
    lgv_run <- tryCatch(require_accepted_upstream("02_local_genetic_variance", a, r),
                        error = function(e) NULL)
    elig <- NA_integer_; excl <- NA_integer_; reasons <- NA_character_
    if (!is.null(lgv_run)) {
        lf <- file.path(V2_ROOT, "02_local_genetic_variance", "_m", "runs",
                        lgv_run$run_id, "results", "combined",
                        sprintf("local-genetic-control-%s-%s-vmrs.tsv", a, r))
        if (file.exists(lf)) {
            l <- fread(lf)
            elig <- sum(l$local_genetic_control_eligible %in% c(TRUE, "TRUE"))
            excl <- nrow(l) - elig
            rr <- l[!(local_genetic_control_eligible %in% c(TRUE, "TRUE")),
                    .N, by = local_genetic_control_exclusion_reason]
            reasons <- paste(sprintf("%s=%d", rr[[1]], rr$N), collapse = "; ")
        }
    }
    data.table(arm = a, region = r,
               vmr_set_id = cat_run$vmr_set_id,
               catalog_run = cat_run$run_id,
               score_run = if (is.null(lgv_run)) NA_character_ else lgv_run$run_id,
               design_n = cohorts$donor_counts[[a]][[r]]$design_n,
               n_donors = nrow(don), n_vmrs_called = nrow(vmr),
               n_vmrs_scored = if (is.na(elig)) NA_integer_ else elig + excl,
               n_vmrs_eligible = elig, n_vmrs_excluded = excl,
               exclusion_reasons = reasons)
}))), fill = TRUE)
write_atomic(denom, file.path(table_dir, "exclusions-and-denominators.tsv"))
message("[table] exclusions-and-denominators.tsv (", nrow(denom), " arm x region cells)")

## ====================================== 4. supplementary table index
##
## The tracked _m/combined/ deliverables a journal would receive.
## `-UNACCEPTED` is excluded BY PATTERN, not by hand: those are outputs of a
## stage whose config is not PI-locked and may not be cited (AGENTS.md 6).
comb <- rbindlist(lapply(MODULES, function(m) {
    d <- file.path(V2_ROOT, m, "_m", "combined")
    if (!dir.exists(d)) return(NULL)
    f <- list.files(d, pattern = "\\.tsv$", full.names = TRUE, recursive = FALSE)
    f <- f[!grepl("-UNACCEPTED|-UNLOCKED", basename(f))]
    if (length(f) == 0) return(NULL)
    data.table(module = m, file = basename(f),
               path = file.path(m, "_m", "combined", basename(f)),
               bytes = file.size(f),
               n_rows = vapply(f, function(x)
                   tryCatch(nrow(fread(x, nrows = Inf)), error = function(e) NA_integer_),
                   integer(1)))
}), fill = TRUE)
write_atomic(comb, file.path(table_dir, "supplementary-table-index.tsv"))
message("[table] supplementary-table-index.tsv (", nrow(comb), " tracked tables)")

## ==================================== 5. software and run manifest
si <- sessioninfo::session_info()
pkgs <- as.data.table(si$packages)[, .(package, loadedversion, source)]
upstream <- sort(unique(unlist(strsplit(
    registry$source_run_id[!is.na(registry$source_run_id)], ";"))))

soft <- rbindlist(list(
    data.table(item = "figure_run_id",  value = opts$run_id),
    data.table(item = "git_commit",     value = git_commit(V2_ROOT)),
    data.table(item = "git_dirty",      value = as.character(git_dirty(V2_ROOT))),
    data.table(item = "r_version",      value = paste(R.version$major, R.version$minor, sep = ".")),
    data.table(item = "conda_prefix",   value = Sys.getenv("CONDA_PREFIX", NA_character_)),
    data.table(item = "hostname",       value = Sys.info()[["nodename"]]),
    data.table(item = "n_figures",      value = as.character(length(rendered))),
    data.table(item = "n_source_tables", value = as.character(nrow(registry))),
    data.table(item = "upstream_runs",  value = paste(upstream, collapse = ";")),
    pkgs[, .(item = paste0("package:", package),
             value = paste0(loadedversion, " [", source, "]"))]))
write_atomic(soft, file.path(table_dir, "software-and-run-manifest.tsv"))
message("[table] software-and-run-manifest.tsv (", length(upstream), " upstream runs)")

## ==================================== 6. booktabs for the claim matrix
##
## Same escaping and structure as 04_table1_cohort.R, generated FROM the TSV so
## the two cannot drift.
esc <- function(x) {
    x <- as.character(x); x[is.na(x)] <- "--"
    gsub("_", "\\\\_", gsub("%", "\\\\%", gsub("&", "\\\\&", x)))
}
tex_cols <- c("module", "run_id", "region", "decision")
tt <- claims[, ..tex_cols]
tex <- c("\\begin{tabular}{llll}", "\\toprule",
         paste0(paste(c("Module", "Accepted run", "Region", "Decision"),
                      collapse = " & "), " \\\\"),
         "\\midrule",
         paste0(apply(tt, 1, function(r) paste(esc(r), collapse = " & ")), " \\\\"),
         "\\bottomrule", "\\end{tabular}", "")
writeLines(tex, file.path(table_dir, "analysis-to-claim-matrix.tex"))
message("[table] analysis-to-claim-matrix.tex")

message("[done] manuscript tables written to ", table_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
