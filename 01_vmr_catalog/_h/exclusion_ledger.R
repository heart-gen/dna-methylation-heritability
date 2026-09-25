#### 01_vmr_catalog / exclusion_ledger: one schema for every exclusion ####
##
## AGENTS.md 7.1 lists a "technical QC and exclusion table" among this module's
## required outputs, 11 says "always report denominators, exclusions, brain
## region, donor group, VMR set", and 14 makes "concordance denominators are
## degenerate or undocumented" a stop condition. Until this file existed the
## module decided exclusions in four places and recorded them in none: donor
## selection and CpG QC in 00_prepare.R, genotype-PC availability in
## 01_analyze.R, the minimum-CpG rule in 02_summarize.R, and the chromosome
## policy in config/thresholds.yml. Every count lived only in a SLURM log.
##
## This file adds no filter and removes none. It is a recording schema, so that
## `qc/exclusions.tsv` can state, for every donor and every candidate VMR that
## entered and did not survive, the single rule that removed it -- and so that
## `qc/exclusion_accounting.tsv` can show entered = survived + excluded rather
## than assert it.
##
## Sourced by 00_prepare.R, 01_analyze.R, 02_summarize.R and 04_turnover.R. Not
## a numbered stage: it emits nothing on its own.

suppressPackageStartupMessages(library(data.table))

## Column order of qc/exclusions.tsv, and of every per-stage part written into a
## run for 04_turnover.R to assemble.
##
##   stage             which script made the decision
##   unit_type         donor | cpg | vmr_candidate | chromosome | phenotype_row
##   unit_id           brnum, "chr1:100-200", "chrX", or NA for aggregate rows
##   chrom             chromosome the decision applied to, or a scope token
##   exclusion_reason  the rule, named so it can be grepped
##   n_units           1 for an itemized row, the count for an aggregate row
##   itemized          FALSE when the units are too numerous to list
EXCLUSION_COLS <- c("stage", "unit_type", "unit_id", "chrom",
                    "exclusion_reason", "n_units", "itemized")

ACCOUNTING_COLS <- c("unit_type", "scope", "n_entered", "n_survived",
                     "n_excluded", "balanced", "note")

#' Build exclusion-ledger rows with the schema enforced.
#'
#' @param stage,unit_type,exclusion_reason character, recycled
#' @param unit_id itemized identifier, or NA for an aggregate row
#' @param n_units 1 for itemized rows; the count for aggregate rows
#' @param itemized FALSE when `unit_id` is NA because the units are not listable
excl_rows <- function(stage, unit_type, exclusion_reason, unit_id = NA_character_,
                      chrom = NA_character_, n_units = 1L, itemized = TRUE) {
    dt <- data.table(
        stage = as.character(stage),
        unit_type = as.character(unit_type),
        unit_id = as.character(unit_id),
        chrom = as.character(chrom),
        exclusion_reason = as.character(exclusion_reason),
        n_units = as.integer(n_units),
        itemized = as.logical(itemized))
    if (nrow(dt) > 0) {
        if (any(is.na(dt$n_units)) || any(dt$n_units < 0)) {
            stop("excl_rows(): n_units must be a non-negative integer")
        }
        if (any(dt$itemized & is.na(dt$unit_id))) {
            stop("excl_rows(): an itemized row needs a unit_id")
        }
        if (any(dt$itemized & dt$n_units != 1L)) {
            stop("excl_rows(): an itemized row counts exactly one unit")
        }
        if (any(!nzchar(dt$exclusion_reason) | is.na(dt$exclusion_reason))) {
            stop("excl_rows(): every excluded unit needs a reason")
        }
    }
    setcolorder(dt, EXCLUSION_COLS)[]
}

#' An empty ledger with the right columns, so rbindlist() is always safe.
excl_empty <- function() excl_rows(character(), character(), character(),
                                   character(), character(), integer(),
                                   logical())

#' First matching reason per element, so each excluded unit has exactly one.
#'
#' @param conditions named list of logical vectors, all the same length, in
#'   precedence order. Names are the reasons.
#' @return character vector; NA where no condition matched (the unit survives).
#'
#' Precedence is what makes the counts add up: a donor who is both under the age
#' floor and absent from the genotype file must be counted once, under the rule
#' that comes first.
first_reason <- function(conditions) {
    if (length(conditions) == 0) stop("first_reason(): no conditions given")
    if (is.null(names(conditions)) || any(!nzchar(names(conditions)))) {
        stop("first_reason(): every condition must be named with its reason")
    }
    n <- unique(vapply(conditions, length, integer(1)))
    if (length(n) != 1L) {
        stop("first_reason(): conditions have differing lengths (",
             paste(vapply(conditions, length, integer(1)), collapse = ", "), ")")
    }
    out <- rep(NA_character_, n)
    for (nm in names(conditions)) {
        hit <- is.na(out) & !is.na(conditions[[nm]]) & conditions[[nm]]
        out[hit] <- nm
    }
    out
}

#' Build one accounting row, and refuse to call an unbalanced ledger balanced.
#'
#' `n_survived` of NA means the survivor count could not be read (for example a
#' QC-refresh run whose source predates the ledger). That is recorded as
#' unbalanced with a note, never silently as balanced.
accounting_row <- function(unit_type, scope, n_entered, n_survived, n_excluded,
                           note = NA_character_) {
    balanced <- !is.na(n_entered) && !is.na(n_survived) && !is.na(n_excluded) &&
        n_entered == n_survived + n_excluded
    data.table(unit_type = unit_type, scope = scope,
               n_entered = as.integer(n_entered),
               n_survived = as.integer(n_survived),
               n_excluded = as.integer(n_excluded),
               balanced = balanced, note = as.character(note))
}

#' Collapse donor rows that are identical on every chromosome into one row.
#'
#' Donor eligibility is a property of the phenotype table, the genotype .psam
#' and the WGBS object, so the same donors are excluded for the same reason on
#' every chromosome. Writing 22 identical rows per donor makes the table
#' unreadable; collapsing only when they really are identical keeps it honest,
#' and a per-chromosome difference stays visible as separate rows.
#' @param by_extra columns beyond the core schema that are part of the identity
#'   of a row. Any column not in the key and not `chrom` must be constant within
#'   a key, or the collapse would silently pick one value; that is checked.
collapse_uniform_chrom <- function(ledger, chroms, scope_token = "autosomes_all",
                                   by_extra = character()) {
    if (nrow(ledger) == 0) return(ledger)
    key_cols <- c(setdiff(EXCLUSION_COLS, "chrom"), by_extra)
    spare <- setdiff(names(ledger), c(key_cols, "chrom"))
    if (length(spare) > 0) {
        stop("collapse_uniform_chrom(): column(s) ", paste(spare, collapse = ", "),
             " are neither part of the key nor the chromosome. Name them in ",
             "by_extra so the collapse cannot drop a value.")
    }
    seen <- ledger[, .(chroms_seen = paste(sort(unique(chrom)), collapse = ",")),
                   by = key_cols]
    all_token <- paste(sort(unique(as.character(chroms))), collapse = ",")
    uniform <- seen[chroms_seen == all_token]
    if (nrow(uniform) == 0) return(ledger)
    keep <- ledger[!uniform, on = key_cols]
    collapsed <- copy(uniform)[, chroms_seen := NULL][, chrom := scope_token]
    setcolorder(collapsed, names(ledger))
    rbindlist(list(collapsed, keep), use.names = TRUE)
}
