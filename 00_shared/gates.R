#### Dependency gates for downstream modules (AGENTS.md 6) ####
##
## "No downstream production run may consume an upstream result until the
## upstream README records a passing acceptance gate and immutable run ID."
##
## Modules 03, 04 and 05 all consume 02_local_genetic_variance, and 05 also
## consumes 01_vmr_catalog directly. Rather than each module re-deriving what
## "accepted" means, the check lives here once.
##
## The gate is deliberately fail-closed and deliberately annoying to bypass. A
## run started with allow_unaccepted = TRUE is a smoke test and is not citable
## (AGENTS.md 14).

#' Parse the "Accepted runs" table out of a module README.
#'
#' The README is the record of record for acceptance (AGENTS.md 6), not a
#' machine-written status file, precisely so that a human has to have looked at
#' the run before anything downstream can consume it. The table is markdown:
#'
#'     | run_id | cohort | region | vmr_set_id | accepted_on | accepted_by |
#'
#' A README whose table is still `_(none)_` yields zero rows, which is the
#' correct answer -- nothing has been accepted.
read_accepted_runs <- function(module_root, root = repo_root()) {
    readme <- file.path(root, module_root, "README.md")
    if (!file.exists(readme)) {
        stop("No README for module '", module_root, "': ", readme)
    }
    lines <- readLines(readme, warn = FALSE)

    ## Take the table under the "Accepted runs" heading, and only that one --
    ## these READMEs contain several other markdown tables.
    start <- grep("^#+\\s*Accepted runs", lines, ignore.case = TRUE)
    if (length(start) == 0) {
        stop("README for '", module_root, "' has no 'Accepted runs' section. ",
             "AGENTS.md 6 requires the acceptance gate be recorded there.")
    }
    tail_lines <- lines[(start[1] + 1):length(lines)]
    nxt <- grep("^#+\\s", tail_lines)
    if (length(nxt) > 0) tail_lines <- tail_lines[seq_len(nxt[1] - 1)]

    rows <- grep("^\\|", tail_lines, value = TRUE)
    ## Drop the header row and the |---|---| separator.
    rows <- rows[!grepl("^\\|[\\s:-]*\\|[\\s:|-]*$", rows, perl = TRUE)]
    if (length(rows) < 2) return(data.table::data.table())

    split_row <- function(r) {
        cells <- strsplit(sub("^\\|", "", sub("\\|\\s*$", "", r)), "|", fixed = TRUE)[[1]]
        trimws(cells)
    }
    header <- split_row(rows[1])
    body <- rows[-1]
    if (length(body) == 0) return(data.table::data.table())

    parsed <- lapply(body, split_row)
    parsed <- parsed[vapply(parsed, length, integer(1)) == length(header)]
    if (length(parsed) == 0) return(data.table::data.table())

    dt <- data.table::as.data.table(
        do.call(rbind, lapply(parsed, function(p) as.list(stats::setNames(p, header)))))
    dt <- dt[!grepl("^_?\\(?none\\)?_?$", dt[[1]], ignore.case = TRUE)]
    dt[]
}

#' Require that an upstream module has accepted a run for this cohort x region.
#'
#' Returns the accepted row (run_id, vmr_set_id, ...) so the caller can record
#' it in its own manifest -- the point of the gate is not only to refuse, but to
#' make the downstream manifest carry the exact upstream handle it consumed
#' (AGENTS.md 9).
#'
#' @param module_root e.g. "02_local_genetic_variance"
#' @param allow_unaccepted TRUE for smoke runs; warns instead of stopping and
#'   returns a row with run_id NA, which callers must propagate into the
#'   manifest as `smoke_run = TRUE`.
require_accepted_upstream <- function(module_root, cohort, region,
                                      allow_unaccepted = FALSE,
                                      root = repo_root()) {
    accepted <- read_accepted_runs(module_root, root = root)

    hit <- if (nrow(accepted) > 0 &&
               all(c("cohort", "region") %in% names(accepted))) {
        accepted[accepted$cohort == cohort & accepted$region == region, ]
    } else {
        accepted[0, ]
    }

    if (nrow(hit) == 1) return(as.list(hit))

    msg <- if (nrow(hit) > 1) {
        paste0("README for '", module_root, "' lists ", nrow(hit),
               " accepted runs for ", cohort, " x ", region,
               ". Exactly one run may be accepted per cell; retire the others.")
    } else {
        paste0("Upstream module '", module_root, "' has no accepted run for ",
               cohort, " x ", region, " (AGENTS.md 6). Record the run ID and ",
               "its passing acceptance gate in ", module_root, "/README.md ",
               "before consuming it downstream.")
    }

    if (allow_unaccepted && nrow(hit) <= 1) {
        warning(msg, "\n  Continuing because allow_unaccepted = TRUE. This run ",
                "is a SMOKE TEST and must not be cited or used in the ",
                "manuscript (AGENTS.md 14).", call. = FALSE)
        return(list(run_id = NA_character_, cohort = cohort, region = region,
                    vmr_set_id = NA_character_, smoke_run = TRUE))
    }
    stop(msg, call. = FALSE)
}

#' Load the manuscript-facing relative local-genetic-control score.
#'
#' The final joint estimator failed its absolute-PVE gate but passed the locked
#' ordering gate. Downstream biological modules therefore consume only the
#' within-cell rank score emitted by Module 02. The raw estimate may remain in
#' the source table for audit/descriptive distributions, but this loader fails
#' unless every row explicitly prohibits absolute-PVE interpretation.
load_local_genetic_control <- function(upstream_run_id, region, cohort,
                                      module_root = "02_local_genetic_variance",
                                      root = repo_root(), eligible_only = TRUE) {
    run_dir <- file.path(root, module_root, "_m", "runs", upstream_run_id)
    if (!dir.exists(run_dir)) {
        stop("Accepted upstream run directory not found: ", run_dir)
    }
    f <- list.files(file.path(run_dir, "results", "combined"),
                    pattern = "^local-genetic-control-.*-vmrs\\.tsv$",
                    full.names = TRUE)
    if (length(f) != 1L) {
        stop("Expected exactly one local-genetic-control table in ", run_dir,
             "/results/combined, found ", length(f))
    }
    dt <- data.table::fread(f[[1L]])
    banned <- intersect(c("h2_unscaled", "r_squared_cv", "h2_en_calibrated"),
                        names(dt))
    if (length(banned)) {
        stop("Relative-score table carries superseded estimator column(s): ",
             paste(banned, collapse = ", "))
    }
    required <- c(
        "vmr_id", "cohort", "region", "vmr_set_id",
        "chrom", "start", "end", "n_cpgs", "n_variants",
        "mean_methylation", "methylation_variance",
        "local_genetic_control_eligible",
        "local_genetic_control_exclusion_reason",
        "local_snp_contribution_score",
        "local_snp_contribution_score_z",
        "local_snp_contribution_quartile",
        "absolute_pve_interpretation_allowed",
        "local_genetic_control_decision"
    )
    missing <- setdiff(required, names(dt))
    if (length(missing)) {
        stop("Relative-score table is missing: ", paste(missing, collapse = ", "))
    }
    if (any(dt$cohort != cohort) || any(dt$region != region)) {
        stop("Relative-score table cohort/region does not match requested cell")
    }
    parse_flag <- function(x, field) {
        value <- tolower(trimws(as.character(x)))
        if (any(is.na(x) | !value %in% c("true", "false", "t", "f", "1", "0"))) {
            stop("Relative-score table has invalid or missing ", field, " values")
        }
        value %in% c("true", "t", "1")
    }
    absolute_allowed <- parse_flag(
        dt$absolute_pve_interpretation_allowed,
        "absolute_pve_interpretation_allowed"
    )
    if (any(absolute_allowed)) {
        stop("Relative-score table authorizes absolute-PVE interpretation")
    }
    if (any(dt$local_genetic_control_decision !=
            "PASS_RELATIVE_GENETIC_CONTROL_FAIL_ABSOLUTE_LOCUS_PVE")) {
        stop("Relative-score table has the wrong interpretation decision")
    }
    eligible <- parse_flag(
        dt$local_genetic_control_eligible,
        "local_genetic_control_eligible"
    )
    if (eligible_only) {
        dt <- dt[eligible]
        eligible <- rep(TRUE, nrow(dt))
    }
    score <- dt$local_snp_contribution_score[eligible]
    if (any(!is.finite(score)) || any(score <= 0 | score >= 1)) {
        stop("Eligible local SNP contribution scores must lie strictly in (0,1)")
    }
    dt[]
}
