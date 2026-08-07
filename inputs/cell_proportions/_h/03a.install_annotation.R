#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(yaml))
ROOT <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
cfg <- yaml::read_yaml(file.path(ROOT, "config/cell_deconvolution.yml"))
ref_dir <- file.path(ROOT, cfg$paths$reference_dir)
lib <- file.path(ref_dir, "R_libs")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

package <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
archive <- file.path(ref_dir, paste0(package, "_", cfg$reference$epic_annotation_version, ".tar.gz"))
if (!file.exists(archive)) stop("Missing checksum-verified annotation archive: ", archive)
if (!requireNamespace(package, quietly = TRUE)) {
  install.packages(archive, repos = NULL, type = "source", lib = lib)
}
stopifnot(as.character(utils::packageVersion(package)) == cfg$reference$epic_annotation_version)
write.table(
  data.frame(package = package, version = as.character(utils::packageVersion(package)),
             library = find.package(package)),
  file.path(ref_dir, "annotation_software_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
