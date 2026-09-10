#!/usr/bin/env Rscript
#### 03 / stage 07: stack the donor-group cells' out-of-fold predictions ####
##
## The two donor-group cells hold DISJOINT donors evaluated on ONE shared VMR
## set, so their per-donor OOF predictions concatenate into a single donor x VMR
## table. That is the "recombine per-individual estimates" half of the design:
## downstream analyses that want every individual back get one table with a
## donor_group column, rather than two files and a positional assumption.
##
## Stacking predictions is well defined here because both cells predict the same
## phenotype (the pooled catalog's VMR mean methylation, from byte-identical
## phenotype files) on the same loci.
##
## Stacking ACCURACY is not. `r2_pred_oof` stays per cell and no pooled r2 is
## computed: 1 - SSE/SST over the union would be inflated by any mean difference
## between the groups, which is a between-group contrast wearing an accuracy
## label. The per-cell values sit side by side instead.
##
## Usage:
##   Rscript _h/07_stack_donor_group_predictions.R --region dlpfc \
##     --run-ids lsp-all_individuals.AA-dlpfc-20260910,lsp-all_individuals.EA-dlpfc-20260910

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

opts <- parse_v2_args(require = c("region", "run_ids"))
region <- tolower(opts$region)
validate_cohort_region(region = region)

run_ids <- trimws(strsplit(opts$run_ids, ",", fixed = TRUE)[[1]])
run_ids <- run_ids[nzchar(run_ids)]
if (length(run_ids) != 2L) {
    stop("Give exactly two run IDs, one per donor group: --run-ids a,b")
}

module_root <- file.path(V2_ROOT, "03_local_snp_prediction")

read_cell <- function(run_id) {
    run_dir <- file.path(module_root, "_m", "runs", run_id)
    if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
    manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
    mval <- function(f) {
        v <- manifest$value[manifest$field == f]
        if (length(v) != 1L) stop("Run ", run_id, " manifest lacks ", f)
        as.character(v[[1L]])
    }
    cell <- mval("cohort")
    parsed <- parse_cell(cell)
    if (!parsed$is_estimation_cell) {
        stop("Run ", run_id, " is discovery arm '", cell, "'. Stacking is for ",
             "the two donor-group cells of one pooled catalog.")
    }
    if (!identical(tolower(mval("region")), region)) {
        stop("Run ", run_id, " is region '", mval("region"), "', not ", region)
    }
    if (!identical(toupper(mval("smoke_run")), "FALSE")) {
        stop("Run ", run_id, " is a smoke run and is not citable")
    }
    accepted <- read_accepted_runs("03_local_snp_prediction")
    if (!nrow(accepted) || !run_id %in% accepted$run_id) {
        stop("Run ", run_id, " is not recorded in the 03 README's accepted ",
             "runs table (AGENTS.md 6)")
    }

    comb <- file.path(run_dir, "results", "combined")
    per_donor <- fread(file.path(comb, "predictions-per-donor.tsv"))
    metrics <- fread(file.path(comb, paste0("oof-prediction-", cell, "-",
                                            region, "-vmrs.tsv")))
    list(run_id = run_id, cell = cell, group = parsed$estimation_group,
         catalog_cohort = parsed$catalog_cohort,
         vmr_set_id = mval("vmr_set_id"),
         per_donor = per_donor, metrics = metrics)
}

cells <- lapply(run_ids, read_cell)
groups <- vapply(cells, function(x) x$group, character(1))
if (anyDuplicated(groups)) {
    stop("Both runs predict in the same donor group (", groups[[1]],
         "). That is not a contrast.")
}
if (length(unique(vapply(cells, function(x) x$catalog_cohort, character(1)))) != 1L) {
    stop("The two cells were discovered on different catalogs")
}
if (length(unique(vapply(cells, function(x) x$vmr_set_id, character(1)))) != 1L) {
    stop("The two cells carry different vmr_set_ids. Discovery must happen ",
         "once in the pooled sample (AGENTS.md 7.7); refusing to stack.")
}

## ------------------------------------------------------- stack the donors
stack_one <- function(x) {
    dt <- copy(as.data.table(x$per_donor))
    ## The column is written by 03_combine_oof.R, but a run predating that
    ## change would have no donor_group; fill it from the manifest rather than
    ## letting the stacked table carry an unlabelled block of rows.
    if (!"donor_group" %in% names(dt)) dt[, donor_group := x$group]
    dt[, `:=`(donor_group = x$group, cohort = x$cell, region = region,
              upstream_lsp_run_id = x$run_id)]
    dt
}
stacked <- rbindlist(lapply(cells, stack_one), use.names = TRUE, fill = TRUE)

## Disjointness is the property that makes a concatenation legitimate. A donor
## in both cells would be double-counted, and would also mean the race filters
## overlap, which 01b's partition proof should already have refused.
dups <- stacked$donor[duplicated(stacked$donor)]
if (length(dups)) {
    stop("Donors appear in both cells: ",
         paste(head(unique(dups), 10), collapse = ", "),
         ". The donor groups must be disjoint.")
}

out_dir <- file.path(module_root, "_m", "combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_atomic(stacked, file.path(
    out_dir, paste0("oof-predictions-per-donor-donor-group-", region, ".tsv")))

## ------------------------------------------- per-cell accuracy, side by side
metric_cols <- c("vmr_id", "n_donors_predicted", "n_predictions",
                 "r2_pred_oof", "cor2_oof", "rmse", "calibration_slope",
                 "screening_pass_frequency", "median_n_variants")
side <- Reduce(function(x, y) merge(x, y, by = "vmr_id", all = FALSE),
               lapply(cells, function(x) {
                   m <- as.data.table(x$metrics)
                   miss <- setdiff(metric_cols, names(m))
                   if (length(miss)) {
                       stop("Run ", x$run_id, " metrics lack: ",
                            paste(miss, collapse = ", "))
                   }
                   m <- m[, metric_cols, with = FALSE]
                   setnames(m, setdiff(metric_cols, "vmr_id"),
                            paste0(setdiff(metric_cols, "vmr_id"), "_", x$group))
                   m
               }))
side[, `:=`(region = region, donor_groups = paste(groups, collapse = ","),
            pooled_r2_emitted = FALSE)]
write_atomic(side, file.path(
    out_dir, paste0("oof-prediction-donor-group-", region, ".tsv")))

summary_dt <- data.table(
    region = region,
    group_a = groups[[1]], group_b = groups[[2]],
    n_donors_a = nrow(cells[[1]]$per_donor),
    n_donors_b = nrow(cells[[2]]$per_donor),
    n_donors_stacked = nrow(stacked),
    n_loci_shared = nrow(side),
    median_r2_a = stats::median(side[[paste0("r2_pred_oof_", groups[[1]])]],
                                na.rm = TRUE),
    median_r2_b = stats::median(side[[paste0("r2_pred_oof_", groups[[2]])]],
                                na.rm = TRUE),
    spearman_r2 = if (nrow(side) > 2) stats::cor(
        side[[paste0("r2_pred_oof_", groups[[1]])]],
        side[[paste0("r2_pred_oof_", groups[[2]])]],
        method = "spearman", use = "complete.obs") else NA_real_,
    pooled_r2_emitted = FALSE
)
print(summary_dt)
write_atomic(summary_dt, file.path(
    out_dir, paste0("donor-group-prediction-summary-", region, ".tsv")))

message("[stack] ", nrow(stacked), " donors across ", nrow(side), " shared loci")
message("  Per-cell r2 only. A pooled r2 over the union would absorb any mean ",
        "difference between the groups; the cells' n and SNP counts differ, so ",
        "compare ordering, not levels (AGENTS.md 7.7).")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
