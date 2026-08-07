#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

ROOT <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
source(file.path(ROOT, "inputs/cell_proportions/_h/dnam_deconvolution_utils.R"))
args <- commandArgs(trailingOnly = TRUE)
region <- normalize_region(if (length(args)) args[[1]] else "caudate")
requested_platform <- if (length(args) >= 2L) toupper(args[[2]]) else "850K"
use_coordinate_wgbs <- requested_platform == "WGBS"
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
cfg <- yaml::read_yaml(file.path(ROOT, "config/cell_deconvolution.yml"))
render <- function(template, region) gsub("\\{region\\}", region, template)
work_template <- if (use_coordinate_wgbs) cfg$paths$dnam_wgbs_work_template else cfg$paths$dnam_work_template
work_dir <- file.path(ROOT, render(work_template, region))
ref_dir <- file.path(ROOT, cfg$paths$reference_dir)
bulk_path <- file.path(work_dir, "scmd_bulk_markers.rds")
if (!file.exists(bulk_path)) stop("Missing extracted bulk markers: ", bulk_path)
bulk <- readRDS(bulk_path)
platform_file <- file.path(ref_dir, if (use_coordinate_wgbs) {
  "scmd_coordinate_wgbs_signature_platform.tsv"
} else {
  "scmd_signature_platform.tsv"
})
signature_platform <- if (file.exists(platform_file)) {
  as.character(fread(platform_file)$signature_platform[[1]])
} else {
  as.character(cfg$reference$primary_signature_platform)
}

local_lib <- file.path(ref_dir, "R_libs")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
set.seed(as.integer(cfg$seed))

engine <- NULL
components <- list()
signature_summary <- data.table()
upstream_error <- NA_character_
result <- NULL

if (requireNamespace("scMD", quietly = TRUE)) {
  result <- tryCatch(
    scMD::scMD(
      bulk = bulk,
      bulk_type = if (use_coordinate_wgbs) "WGBS" else cfg$engine$bulk_type,
      ncluster = n_cores,
      nmrk = as.integer(cfg$reference$markers_per_cell),
      enableFileSaving = FALSE
    ),
    error = function(e) e
  )
  if (!inherits(result, "error")) {
    engine <- "scMD_1.0.0_upstream"
    expected_cells <- c("Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc")
    ensemble <- orient_fraction_matrix(result$scMD_p, colnames(bulk), expected_cells)
    components <- result$phat_all
    saveRDS(result$phat_all, file.path(work_dir, "scmd_upstream_phat_all.rds"))
  } else {
    upstream_error <- conditionMessage(result)
  }
}

if (is.null(engine)) {
  if (!isTRUE(cfg$engine$allow_reference_ensemble_fallback)) {
    stop("Upstream scMD unavailable or failed: ", upstream_error)
  }
  refs <- if (use_coordinate_wgbs) {
    list(
      Lee = load_reference_objects(file.path(ref_dir, "Lee_7ct_WGBS.rda"), "Lee", "WGBS"),
      Tian = load_reference_objects(file.path(ref_dir, "Tian_7ct_WGBS.rda"), "Tian", "WGBS")
    )
  } else {
    list(
      Lee = load_reference_objects(file.path(ref_dir, "Lee_7ct_450850.rda"), "Lee", "450850"),
      Tian = load_reference_objects(file.path(ref_dir, "Tian_7ct_450850.rda"), "Tian", "450850")
    )
  }
  result <- run_reference_ensemble(
    bulk, refs, as.integer(cfg$reference$markers_per_cell)
  )
  ensemble <- result$ensemble
  components <- result$components
  signature_summary <- result$signature_summary
  engine <- cfg$engine$fallback_name
  saveRDS(components, file.path(work_dir, "scmd_reference_components.rds"))
}

labels <- cfg$cell_type_labels
long <- as.data.table(ensemble, keep.rownames = "sample_id")
long <- melt(long, id.vars = "sample_id", variable.name = "cell_type_raw",
             value.name = "proportion")
long[, `:=`(
  cell_type = standardize_cell_type_names(cell_type_raw, labels),
  region = region,
  method = engine,
  reference = paste0("Lee_Tian_", signature_platform)
)]
setcolorder(long, c("sample_id", "cell_type", "proportion", "region", "method", "reference", "cell_type_raw"))

component_rows <- rbindlist(lapply(names(components), function(nm) {
  x <- tryCatch(
    orient_fraction_matrix(components[[nm]], colnames(bulk), colnames(ensemble)),
    error = function(e) NULL
  )
  if (is.null(x)) return(NULL)
  d <- as.data.table(x, keep.rownames = "sample_id")
  d <- melt(d, id.vars = "sample_id", variable.name = "cell_type_raw", value.name = "proportion")
  d[, `:=`(
    cell_type = standardize_cell_type_names(cell_type_raw, labels),
    region = region,
    component = nm
  )]
  d
}), fill = TRUE)

prop_template <- if (use_coordinate_wgbs) cfg$paths$dnam_wgbs_proportions_template else cfg$paths$dnam_proportions_template
out_prop <- file.path(ROOT, render(prop_template, region))
fwrite(long, out_prop, sep = "\t")
fwrite(component_rows, file.path(work_dir, "dnam-scmd-component-proportions.tsv.gz"), sep = "\t")
if (nrow(signature_summary)) {
  fwrite(signature_summary, file.path(work_dir, "signature_summary.tsv"), sep = "\t")
}

row_sums <- rowSums(ensemble)
qc <- data.table(
  region = region,
  sample_id = rownames(ensemble),
  n_cell_types = ncol(ensemble),
  fraction_sum = row_sums,
  min_fraction = apply(ensemble, 1, min),
  max_fraction = apply(ensemble, 1, max),
  bounded = apply(ensemble, 1, function(x) all(is.finite(x) & x >= 0 & x <= 1)),
  sum_to_one = abs(row_sums - 1) <= as.numeric(cfg$validation$fraction_sum_tolerance)
)
fwrite(qc, file.path(work_dir, "deconvolution_qc.tsv"), sep = "\t")
manifest <- data.table(
  region = region,
  engine = engine,
  upstream_scmd_installed = requireNamespace("scMD", quietly = TRUE),
  upstream_error = upstream_error,
  n_samples = nrow(ensemble),
  n_cell_types = ncol(ensemble),
  n_bulk_marker_ids = nrow(bulk),
  signature_platform = signature_platform,
  seed = as.integer(cfg$seed)
)
fwrite(manifest, file.path(work_dir, "deconvolution_engine_manifest.tsv"), sep = "\t", na = "NA")
capture.output(sessionInfo(), file = file.path(work_dir, "deconvolution_session_info.txt"))
print(manifest)
