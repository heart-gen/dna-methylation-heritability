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

    ## rbind() over as.list() rows yields a matrix of lists, and
    ## as.data.table() then keeps every column as a list column, so the
    ## `accepted$cohort == cohort` comparison in require_accepted_upstream()
    ## fails with "comparison of these types is not implemented". Build plain
    ## character columns instead -- every cell here is markdown text.
    dt <- data.table::as.data.table(stats::setNames(
        lapply(seq_along(header),
               function(j) vapply(parsed, `[`, character(1), j)),
        header))
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

    ## Compute the row selection OUTSIDE the data.table `[`. Inside it, the
    ## bare names `cohort` and `region` resolve to the COLUMNS rather than to
    ## these arguments, so `accepted$cohort == cohort` is column == column,
    ## which is TRUE for every row and matches the whole table.
    hit <- if (nrow(accepted) > 0 &&
               all(c("cohort", "region") %in% names(accepted))) {
        keep <- accepted$cohort == cohort & accepted$region == region
        accepted[which(keep), ]
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

#' The donor-group inference policy, read from config rather than hard-coded.
#'
#' `config/analysis_thresholds.yml:donor_group` states what the AA-vs-EA axis is
#' allowed to conclude. It replaced `require_interaction_for_ancestry_claim`,
#' which named a test the design cannot run -- an interaction term needs a
#' pooled group x genotype model and AGENTS.md 7.6 forbids comparing the Module
#' 02 score across cells at any level -- and so guarded nothing.
#'
#' Every stage that touches both cells calls this and asserts against it, so the
#' restriction lives in one place and a config edit cannot quietly widen a claim
#' without a stage refusing. Returns the three flags plus the stratification
#' floor.
donor_group_inference_policy <- function(root = repo_root()) {
    cfg <- load_config("analysis_thresholds", root = root)
    dg <- config_get(cfg, "donor_group")

    required <- c("donor_group_inference", "cross_group_raw_score_comparison",
                  "ancestry_effect_claim_allowed",
                  "min_n_per_group_for_stratified")
    missing <- setdiff(required, names(dg))
    if (length(missing)) {
        stop("config/analysis_thresholds.yml:donor_group lacks: ",
             paste(missing, collapse = ", "),
             ". The donor-group axis will not run against an underspecified ",
             "policy -- state every key explicitly (AGENTS.md 7.7).")
    }

    ## `concordance_only` is the only implemented mode. Anything else is a
    ## scientific decision that needs the eliminating analysis and new code, so
    ## refuse rather than fall through to a default.
    mode <- as.character(dg$donor_group_inference)
    if (!identical(mode, "concordance_only")) {
        stop("Unsupported donor_group_inference '", mode, "'. Only ",
             "'concordance_only' is implemented; the reportable quantity is ",
             "ordering agreement on the shared locus set, never a level ",
             "(AGENTS.md 7.6, 7.7).")
    }

    flag <- function(field) {
        v <- dg[[field]]
        if (!is.logical(v) || length(v) != 1L || is.na(v)) {
            stop("config/analysis_thresholds.yml:donor_group$", field,
                 " must be true or false, got: ", paste(v, collapse = ", "))
        }
        v
    }
    raw_ok <- flag("cross_group_raw_score_comparison")
    ancestry_ok <- flag("ancestry_effect_claim_allowed")

    ## These two are consequences of the locked rank scope and of 7.7's
    ## elimination requirement, not tunables. Catch a well-meaning edit here
    ## rather than in a figure caption.
    if (raw_ok) {
        stop("cross_group_raw_score_comparison is true. The Module 02 score is ",
             "a within-cell midrank percentile and ",
             "config/local_genetic_control.yml locks rank_scope: ",
             "cohort_by_region, so cross-cell levels are not on one scale ",
             "(AGENTS.md 7.6).")
    }
    if (ancestry_ok) {
        stop("ancestry_effect_claim_allowed is true. Sample size, MAF, LD and ",
             "SNP availability differ between the cells and must be eliminated ",
             "before any difference is attributed to ancestry (AGENTS.md 7.7). ",
             "Flipping this key is not the eliminating analysis.")
    }

    floor_n <- suppressWarnings(as.integer(dg$min_n_per_group_for_stratified))
    if (!is.finite(floor_n) || floor_n < 1L) {
        stop("min_n_per_group_for_stratified must be a positive integer")
    }

    list(mode = mode,
         cross_group_raw_score_comparison = raw_ok,
         ancestry_effect_claim_allowed = ancestry_ok,
         min_n_per_group_for_stratified = floor_n,
         config_sha256 = attr(cfg, "config_sha256"))
}
