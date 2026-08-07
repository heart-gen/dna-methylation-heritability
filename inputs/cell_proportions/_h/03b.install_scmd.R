#!/usr/bin/env Rscript

# Install scMD and its two GitHub dependencies into a project-local library.
# All revisions are locked in config/cell_deconvolution.yml. The source archive
# is checksum-verified by 03.download_scmd_reference.sh before this script runs.

suppressPackageStartupMessages({
  library(yaml)
  library(remotes)
})

ROOT <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
cfg <- yaml::read_yaml(file.path(ROOT, "config/cell_deconvolution.yml"))
ref_dir <- file.path(ROOT, cfg$paths$reference_dir)
lib <- file.path(ref_dir, "R_libs")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

# Fail before network installation when the active R build cannot compile even
# a minimal dependency. This was observed in the legacy epigenomics prefix.
config_test <- suppressWarnings(system2(file.path(R.home("bin"), "R"),
  c("CMD", "config", "SHLIB_LIBADD"), stdout = TRUE, stderr = TRUE))
if (!is.null(attr(config_test, "status")) && attr(config_test, "status") != 0L) {
  stop("Active R installation lacks Makeconf variable SHLIB_LIBADD; create the locked cell-deconvolution conda environment before installing upstream scMD")
}

deps <- cfg$reference$direct_dependency_commits
install_locked <- function(package, repo, commit) {
  current <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
  remotes::install_github(
    paste0(repo, "@", commit), lib = lib, dependencies = TRUE,
    upgrade = "never", build_vignettes = FALSE, quiet = FALSE
  )
  data.frame(package = package, repository = repo, commit = commit,
             prior_version = current, installed_version = as.character(utils::packageVersion(package)))
}

records <- list(
  install_locked("EnsDeconv", "randel/EnsDeconv", deps$EnsDeconv),
  install_locked("MIND", "randel/MIND", deps$MIND)
)

archive <- file.path(ref_dir, paste0("scMD-", cfg$reference$scmd_commit, ".tar.gz"))
if (!file.exists(archive)) stop("Missing checksum-verified scMD source archive: ", archive)
remotes::install_local(
  archive, lib = lib, dependencies = TRUE, upgrade = "never",
  build_vignettes = FALSE, quiet = FALSE
)
records[[length(records) + 1L]] <- data.frame(
  package = "scMD", repository = cfg$reference$scmd_tarball_url,
  commit = cfg$reference$scmd_commit, prior_version = NA_character_,
  installed_version = as.character(utils::packageVersion("scMD"))
)

manifest <- do.call(rbind, records)
write.table(manifest, file.path(ref_dir, "scmd_software_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(ref_dir, "scmd_install_session_info.txt"))
print(manifest)
