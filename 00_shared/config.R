#### Configuration loading and path resolution (v2 revision) ####
##
## AGENTS.md 9: "Quest paths belong in configuration or environment variables,
## never inside analysis functions." Every path a v2 script touches is resolved
## through resolve_path() against config/, so the six-way copy-paste drift that
## produced defects V2, V3 and V11 cannot recur.
##
## AGENTS.md 12/14: scientific decisions the PI has not locked must not be made
## silently. assert_locked() stops a production run that would consume an
## unlocked key; smoke runs opt out with allow_unlocked = TRUE.

suppressPackageStartupMessages({
    library(yaml)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Locate the repository root the way here::here() does.
#'
#' Walks up from `start` until it finds a directory containing .git. Used by the
#' shell launchers too (00_shared/slurm.sh), so both agree on the root.
repo_root <- function(start = getwd()) {
    dir <- normalizePath(start, mustWork = TRUE)
    while (dir != dirname(dir)) {
        if (dir.exists(file.path(dir, ".git"))) return(dir)
        dir <- dirname(dir)
    }
    stop("Could not locate repository root (no .git found above ", start, ")")
}

#' Load a YAML config from config/ by bare name.
#'
#' @param name e.g. "paths", "cohorts", "thresholds"
load_config <- function(name, root = repo_root()) {
    f <- file.path(root, "config", paste0(name, ".yml"))
    if (!file.exists(f)) stop("Config not found: ", f)
    cfg <- yaml::read_yaml(f)
    attr(cfg, "config_file") <- f
    attr(cfg, "config_sha256") <- file_sha256(f)
    cfg
}

#' SHA-256 of a file, for the run manifest (AGENTS.md 9).
file_sha256 <- function(path) {
    if (!file.exists(path)) return(NA_character_)
    if (requireNamespace("digest", quietly = TRUE)) {
        return(digest::digest(path, algo = "sha256", file = TRUE))
    }
    ## digest is not in every env; fall back to the system tool rather than
    ## silently recording NA for a provenance field.
    out <- tryCatch(system2("sha256sum", shQuote(path), stdout = TRUE),
                    error = function(e) NA_character_)
    if (length(out) == 0 || is.na(out[1])) return(NA_character_)
    sub(" .*$", "", out[1])
}

#' Fetch a nested config value by dotted key, erroring rather than returning NULL.
config_get <- function(cfg, key) {
    parts <- strsplit(key, ".", fixed = TRUE)[[1]]
    val <- cfg
    for (p in parts) {
        if (!is.list(val) || !p %in% names(val)) {
            stop("Config key not found: ", key)
        }
        val <- val[[p]]
    }
    val
}

#' Resolve a path template from config/paths.yml.
#'
#' Templates use {region}, {chrom}, {cohort}, {run_id}, {module_root}. Relative
#' results are made absolute against the repository root; absolute values in the
#' config (resources on Quest outside the repo) are returned untouched.
#'
#' @param key dotted key into paths.yml, e.g. "wgbs_bsobj_template"
#' @param check if TRUE, stop when the resolved path does not exist
resolve_path <- function(key, region = NULL, chrom = NULL, cohort = NULL,
                         run_id = NULL, module_root = NULL,
                         paths = NULL, root = repo_root(), check = FALSE) {
    if (is.null(paths)) paths <- load_config("paths", root = root)
    tmpl <- config_get(paths, key)
    if (is.null(tmpl) || (length(tmpl) == 1 && is.na(tmpl))) {
        stop("Path key is null in config/paths.yml: ", key)
    }
    subs <- list(region = region, chrom = chrom, cohort = cohort,
                 run_id = run_id, module_root = module_root)
    for (nm in names(subs)) {
        token <- paste0("{", nm, "}")
        if (grepl(token, tmpl, fixed = TRUE)) {
            if (is.null(subs[[nm]])) {
                stop("Path template '", key, "' needs ", nm, " but none was given")
            }
            tmpl <- gsub(token, as.character(subs[[nm]]), tmpl, fixed = TRUE)
        }
    }
    leftover <- regmatches(tmpl, regexpr("\\{[a-z_]+\\}", tmpl))
    if (length(leftover) > 0) {
        stop("Unsubstituted token ", leftover, " in path template '", key, "'")
    }
    out <- if (startsWith(tmpl, "/")) tmpl else file.path(root, tmpl)
    if (check && !file.exists(out)) {
        stop("Resolved path does not exist: ", out, " (key '", key, "')")
    }
    out
}

#' Refuse to start a production run on unlocked PI decisions.
#'
#' AGENTS.md 12: "Agents may recommend defaults but must not silently make these
#' scientific decisions." AGENTS.md 14 makes an unlocked primary cohort,
#' chromosome policy, or VMR covariate set a stop condition.
#'
#' @param cfgs named list of configs the run will consume
#' @param allow_unlocked TRUE for smoke runs; warns instead of stopping
assert_locked <- function(cfgs, allow_unlocked = FALSE) {
    if (!is.list(cfgs) || is.null(names(cfgs))) {
        stop("assert_locked() expects a named list of configs")
    }
    unlocked <- names(cfgs)[!vapply(cfgs, function(c) isTRUE(c$pi_locked), logical(1))]
    if (length(unlocked) == 0) return(invisible(TRUE))

    msg <- paste0(
        "PI decisions are not locked in: ",
        paste(unlocked, collapse = ", "),
        ". Set pi_locked: true in each config once the PI has signed off ",
        "(AGENTS.md 12)."
    )
    if (allow_unlocked) {
        warning(msg, "\n  Continuing because allow_unlocked = TRUE. ",
                "This run is a SMOKE TEST and must not be used as production ",
                "or cited downstream.", call. = FALSE)
        return(invisible(FALSE))
    }
    stop(msg, "\n  Pass --allow-unlocked to run a smoke test.", call. = FALSE)
}

#' Standard command-line parsing for v2 analysis scripts.
#'
#' Every 01_vmr_catalog script takes the same arguments, which is what makes one
#' codepath serve both cohorts and all three regions.
parse_v2_args <- function(args = commandArgs(trailingOnly = TRUE),
                          require = c("cohort", "region")) {
    opts <- list(allow_unlocked = FALSE)
    i <- 1
    while (i <= length(args)) {
        a <- args[[i]]
        if (a == "--allow-unlocked") {
            opts$allow_unlocked <- TRUE
            i <- i + 1
        } else if (startsWith(a, "--")) {
            nm <- gsub("-", "_", sub("^--", "", a))
            if (i + 1 > length(args)) stop("Missing value for ", a)
            opts[[nm]] <- args[[i + 1]]
            i <- i + 2
        } else {
            stop("Unexpected positional argument: ", a,
                 ". v2 scripts take named arguments only ",
                 "(--cohort, --region, --chrom, --run-id).")
        }
    }
    missing <- setdiff(require, names(opts))
    if (length(missing) > 0) {
        stop("Missing required argument(s): --", paste(missing, collapse = " --"))
    }
    validate_cohort_region(opts$cohort, opts$region)
    opts
}

validate_cohort_region <- function(cohort = NULL, region = NULL, root = repo_root()) {
    cohorts <- load_config("cohorts", root = root)
    if (!is.null(cohort) && !cohort %in% cohorts$arms) {
        stop("Unknown cohort '", cohort, "'. Valid: ",
             paste(cohorts$arms, collapse = ", "))
    }
    if (!is.null(region) && !region %in% cohorts$regions) {
        stop("Unknown region '", region, "'. Valid: ",
             paste(cohorts$regions, collapse = ", "))
    }
    invisible(TRUE)
}

#' Definition block for one cohort arm, with paths already made absolute.
cohort_def <- function(cohort, root = repo_root()) {
    cohorts <- load_config("cohorts", root = root)
    validate_cohort_region(cohort = cohort, root = root)
    d <- cohorts$arm_definitions[[cohort]]
    for (k in c("phenotype_table", "psam", "pgen_prefix")) {
        if (!is.null(d[[k]]) && !startsWith(d[[k]], "/")) {
            d[[k]] <- file.path(root, d[[k]])
        }
    }
    d$cohort <- cohort
    d
}

#' Per-region sample blacklist, or NULL.
#'
#' Every region is null in v2: the legacy blacklists only reconciled a stale AA
#' phenotype file and are retired (see config/cohorts.yml). This function is
#' kept so the guard it enforces survives -- it never falls back to another
#' region's file. The legacy hippocampus copy read the DLPFC blacklist
#' (all_individuals/hippocampus/_h/02b.res_var.R:81), silently excluding the
#' wrong donors.
sample_blacklist <- function(region, root = repo_root()) {
    cohorts <- load_config("cohorts", root = root)
    p <- cohorts$sample_blacklist[[region]]
    if (is.null(p)) return(NULL)
    f <- if (startsWith(p, "/")) p else file.path(root, p)
    if (!file.exists(f)) {
        stop("Blacklist declared for region '", region, "' but not found: ", f,
             "\n  Set it to null in config/cohorts.yml if there are no exclusions.")
    }
    ids <- readLines(f, warn = FALSE)
    trimws(ids[nzchar(trimws(ids))])
}
