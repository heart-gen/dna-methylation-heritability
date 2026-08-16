#### Bootstrap for the v2 shared library ####
##
## Every v2 analysis script starts with:
##
##     source(file.path(Sys.getenv("V2_ROOT", "."), "00_shared", "load.R"))
##
## or, when the working directory is anywhere inside the repo:
##
##     source("00_shared/load.R")   # after cd'ing to the root
##
## This file must have no dependencies of its own, because it is what makes the
## rest available. It self-locates the repository root so scripts can be invoked
## from a module's _m/ directory (the submission convention) without every
## script re-deriving paths.

local({
    find_root <- function(start) {
        dir <- normalizePath(start, mustWork = TRUE)
        while (dir != dirname(dir)) {
            if (dir.exists(file.path(dir, ".git"))) return(dir)
            dir <- dirname(dir)
        }
        stop("Could not locate repository root (no .git above ", start, ")")
    }

    root <- Sys.getenv("V2_REPO_ROOT", "")
    if (!nzchar(root)) root <- find_root(getwd())

    shared <- file.path(root, "00_shared")
    if (!dir.exists(shared)) {
        stop("00_shared/ not found under ", root)
    }

    ## Order matters: config.R defines repo_root() and load_config(), which the
    ## others call at load time.
    for (f in c("config.R", "identity.R", "chrom.R", "runid.R", "wgbs.R")) {
        source(file.path(shared, f), local = FALSE)
    }

    assign("V2_ROOT", root, envir = globalenv())
    Sys.setenv(V2_REPO_ROOT = root)
})
