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
    ## --group is the donor-group axis (AGENTS.md 7.7). It composes the cell
    ## token rather than travelling as a second argument, so exactly one token
    ## identifies the cell everywhere downstream. Passing an already-composed
    ## --cohort all_individuals.EA is equivalent and also accepted.
    if (!is.null(opts$group) && nzchar(opts$group)) {
        if (is.null(opts$cohort)) {
            stop("--group requires --cohort (the discovery arm)")
        }
        if (grepl(".", opts$cohort, fixed = TRUE)) {
            stop("Pass either --cohort <arm> --group <group> or a composed ",
                 "--cohort <arm>.<group>, not both: ", opts$cohort,
                 " + ", opts$group)
        }
        opts$cohort <- paste0(opts$cohort, ".", opts$group)
    }
    validate_cohort_region(opts$cohort, opts$region)
    if (!is.null(opts$cohort)) {
        parsed <- parse_cell(opts$cohort)
        opts$cell <- parsed$cell
        opts$catalog_cohort <- parsed$catalog_cohort
        opts$estimation_group <- parsed$estimation_group
    }
    opts
}

#' Every valid cell token: the discovery arms plus the estimation cells.
#'
#' `cohort` in v2 means "the cell this run is identified by". For the arms that
#' is the discovery population; for an estimation cell it is
#' "{catalog_cohort}.{estimation_group}". Both live in the same namespace so the
#' token can go wherever the cohort token went -- run IDs, README acceptance
#' tables, gates.R, downstream manifest fields -- with no schema change.
valid_cells <- function(root = repo_root()) {
    cohorts <- load_config("cohorts", root = root)
    c(as.character(cohorts$arms), names(cohorts$estimation_cells))
}

validate_cohort_region <- function(cohort = NULL, region = NULL, root = repo_root()) {
    cohorts <- load_config("cohorts", root = root)
    cells <- valid_cells(root = root)
    if (!is.null(cohort) && !cohort %in% cells) {
        stop("Unknown cohort/cell '", cohort, "'. Valid: ",
             paste(cells, collapse = ", "))
    }
    if (!is.null(region) && !region %in% cohorts$regions) {
        stop("Unknown region '", region, "'. Valid: ",
             paste(cohorts$regions, collapse = ", "))
    }
    invisible(TRUE)
}

#' Split a cell token into its discovery arm and its estimation group.
#'
#' A BARE arm token parses as a cell whose estimation group equals its cohort --
#' that is what every pre-2026-09-10 run is, and why widening the namespace does
#' not change any existing behaviour.
#'
#' @return list(cell, catalog_cohort, estimation_group, is_estimation_cell)
parse_cell <- function(cell, root = repo_root()) {
    cohorts <- load_config("cohorts", root = root)
    cell <- as.character(cell)
    if (length(cell) != 1L || is.na(cell) || !nzchar(cell)) {
        stop("parse_cell() needs one non-empty cell token")
    }
    if (cell %in% as.character(cohorts$arms)) {
        return(list(cell = cell, catalog_cohort = cell, estimation_group = cell,
                    is_estimation_cell = FALSE))
    }
    spec <- cohorts$estimation_cells[[cell]]
    if (is.null(spec)) {
        stop("Unknown cohort/cell '", cell, "'. Valid: ",
             paste(valid_cells(root = root), collapse = ", "))
    }
    for (k in c("catalog_cohort", "estimation_group", "race_filter")) {
        if (is.null(spec[[k]])) {
            stop("estimation_cells:", cell, " lacks ", k, " in config/cohorts.yml")
        }
    }
    ## The catalog a cell is discovered on must itself be a real arm, or the
    ## cell would point at a Module 01 run that cannot exist.
    if (!spec$catalog_cohort %in% as.character(cohorts$arms)) {
        stop("estimation_cells:", cell, " names catalog_cohort '",
             spec$catalog_cohort, "', which is not an arm")
    }
    ## Guard the token convention itself, so a mislabelled key cannot silently
    ## send Module 02 at the wrong donors.
    expected <- paste0(spec$catalog_cohort, ".", spec$estimation_group)
    if (!identical(cell, expected)) {
        stop("estimation_cells key '", cell, "' does not match its own ",
             "definition; expected '", expected, "'")
    }
    list(cell = cell, catalog_cohort = spec$catalog_cohort,
         estimation_group = spec$estimation_group, is_estimation_cell = TRUE)
}

#' Full definition for one cell: catalog paths from its arm, donors from itself.
#'
#' For a bare arm this returns cohort_def() with the cell fields added, so the
#' two are interchangeable at the call site.
cell_def <- function(cell, root = repo_root()) {
    parsed <- parse_cell(cell, root = root)
    d <- cohort_def(parsed$catalog_cohort, root = root)
    if (parsed$is_estimation_cell) {
        cohorts <- load_config("cohorts", root = root)
        spec <- cohorts$estimation_cells[[cell]]
        ## The cell narrows the donor set. It inherits the catalog arm's
        ## genotype and phenotype FILES -- an estimation cell is a subset of a
        ## pooled run, never a different genotype build.
        d$race_filter <- as.character(spec$race_filter)
        d$label <- spec$label
    }
    d$cohort <- parsed$catalog_cohort
    d$cell <- parsed$cell
    d$catalog_cohort <- parsed$catalog_cohort
    d$estimation_group <- parsed$estimation_group
    d$is_estimation_cell <- parsed$is_estimation_cell
    d
}

#' Definition block for one cohort arm, with paths already made absolute.
cohort_def <- function(cohort, root = repo_root()) {
    cohorts <- load_config("cohorts", root = root)
    validate_cohort_region(cohort = cohort, root = root)
    if (!cohort %in% as.character(cohorts$arms)) {
        stop("cohort_def() takes a discovery ARM, not the estimation cell '",
             cohort, "'. Use cell_def(), or parse_cell()$catalog_cohort.")
    }
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
