#### File checksums (AGENTS.md 9) ####
##
## `file_sha256()` lived in config.R until 2026-09-18. That file does
## `library(yaml)`, and the stages that run in the `calibrated-local-h2`
## estimator env cannot load yaml -- the env is kept byte-for-byte what the
## frozen model was fitted with, which is why `submit_observed_local_control.sh`
## splits Stage 00 (epigenomics) from Stages 01-06 (calibrated-local-h2).
##
## 02/_h/03_apply_frozen_joint_model.R needs exactly this one function, to verify
## that the characterized-support table it scores against is the one Stage 00
## pinned in the manifest. It sources only its module's `00_functions.R`, so the
## call resolved to nothing and the stage died with
## `could not find function "file_sha256"` -- after the 2,306-task feature array
## had already run. Keeping the definition here, with no dependencies, lets both
## envs share ONE implementation; a second copy in the module would be a
## provenance function that could drift from the one writing the manifest.

#' SHA-256 of a file, for the run manifest (AGENTS.md 9).
file_sha256 <- function(path) {
    if (!file.exists(path)) return(NA_character_)
    if (requireNamespace("digest", quietly = TRUE)) {
        return(digest::digest(path, algo = "sha256", file = TRUE))
    }
    ## digest is not in every env; fall back to the system tool rather than
    ## silently recording NA for a provenance field. Both paths hash the file
    ## bytes, so the two envs agree on the value.
    out <- tryCatch(system2("sha256sum", shQuote(path), stdout = TRUE),
                    error = function(e) NA_character_)
    if (length(out) == 0 || is.na(out[1])) return(NA_character_)
    sub(" .*$", "", out[1])
}
