#### Run identity, provenance, and task reconciliation (v2 revision) ####
##
## AGENTS.md 9: every production output must carry run ID, git commit, config
## checksum, input checksums, upstream run IDs, vmr_set_id, ordered donor
## checksum, region/cohort/n, seeds, environment, SLURM job IDs, and output
## checksums. Runs are immutable: "Never update a completed run in place."
##
## AGENTS.md 9 also: "Production runs have zero tolerance for unexplained
## computational failures." reconcile() is what enforces that -- a SLURM array
## that exits 0 on every task is not evidence that every task produced output.

suppressPackageStartupMessages({
    library(data.table)
})

#' Build a run ID: {module}-{cohort}-{region}-{YYYYMMDD}, with a letter suffix
#' if that directory already exists.
#'
#' Never reuses or overwrites an existing run directory.
make_run_id <- function(module, cohort, region, module_root, date = Sys.Date()) {
    base <- sprintf("%s-%s-%s-%s", module, cohort, region, format(date, "%Y%m%d"))
    runs_dir <- file.path(module_root, "_m", "runs")
    candidate <- base
    suffix <- letters
    i <- 1
    while (dir.exists(file.path(runs_dir, candidate))) {
        if (i > length(suffix)) {
            stop("More than ", length(suffix), " runs for ", base,
                 " in one day; something is looping.")
        }
        candidate <- paste0(base, "-", suffix[[i]])
        i <- i + 1
    }
    candidate
}

#' Create an immutable run directory and write its provenance manifest.
#'
#' @param upstream named list of upstream run IDs, e.g. list(vmr_catalog = "...")
#' @return list with run_id, dir, and the manifest fields
new_run <- function(module, cohort, region, module_root,
                    run_id = NULL, upstream = list(), vmr_set_id = NA_character_,
                    extra = list(), root = repo_root()) {
    if (is.null(run_id)) {
        run_id <- make_run_id(module, cohort, region, module_root)
    }
    run_dir <- file.path(module_root, "_m", "runs", run_id)
    if (dir.exists(run_dir)) {
        stop("Run directory already exists and runs are immutable: ", run_dir,
             "\n  AGENTS.md 5.2: never update a completed run in place.")
    }
    dir.create(file.path(run_dir, "logs"), recursive = TRUE)
    dir.create(file.path(run_dir, "excluded"), recursive = TRUE, showWarnings = FALSE)

    manifest <- c(
        list(
            run_id          = run_id,
            analysis        = module,
            cohort          = cohort,
            region          = region,
            vmr_set_id      = vmr_set_id,
            git_commit      = git_commit(root),
            git_dirty       = git_dirty(root),
            started_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
            r_version       = paste(R.version$major, R.version$minor, sep = "."),
            conda_prefix    = Sys.getenv("CONDA_PREFIX", NA_character_),
            hostname        = Sys.info()[["nodename"]],
            slurm_job_id    = Sys.getenv("SLURM_JOB_ID", NA_character_),
            slurm_array_id  = Sys.getenv("SLURM_ARRAY_JOB_ID", NA_character_),
            config_paths_sha256      = file_sha256(file.path(root, "config", "paths.yml")),
            config_cohorts_sha256    = file_sha256(file.path(root, "config", "cohorts.yml")),
            config_thresholds_sha256 = file_sha256(file.path(root, "config", "thresholds.yml")),
            config_covariates_sha256 = file_sha256(file.path(root, "config", "covariates.yml"))
        ),
        ## Guarded: paste0("upstream_", NULL) returns "upstream_", a length-1
        ## vector, so an empty upstream list would produce a name with no value.
        if (length(upstream) > 0) {
            stats::setNames(as.list(unlist(upstream)),
                            paste0("upstream_", names(upstream)))
        } else list(),
        extra
    )

    write_manifest(run_dir, manifest)
    message("[run] created ", run_dir)
    list(run_id = run_id, dir = run_dir, module_root = module_root,
         manifest = manifest)
}

write_manifest <- function(run_dir, manifest) {
    dt <- data.table::data.table(
        field = names(manifest),
        value = vapply(manifest, function(v) {
            if (is.null(v) || length(v) == 0) NA_character_ else as.character(v)[1]
        }, character(1))
    )
    write_atomic(dt, file.path(run_dir, "manifest.tsv"))
}

#' Append fields to a run manifest that is still open (before the run closes).
append_manifest <- function(run, fields) {
    f <- file.path(run$dir, "manifest.tsv")
    dt <- data.table::fread(f, colClasses = "character")
    add <- data.table::data.table(
        field = names(fields),
        value = vapply(fields, function(v) {
            if (is.null(v) || length(v) == 0) NA_character_ else as.character(v)[1]
        }, character(1))
    )
    dt <- rbind(dt[!field %in% add$field], add)
    write_atomic(dt, f)
    invisible(dt)
}

git_commit <- function(root = repo_root()) {
    out <- tryCatch(
        system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"),
                stdout = TRUE, stderr = FALSE),
        error = function(e) NA_character_)
    if (length(out) == 0) NA_character_ else out[1]
}

git_dirty <- function(root = repo_root()) {
    out <- tryCatch(
        system2("git", c("-C", shQuote(root), "status", "--porcelain"),
                stdout = TRUE, stderr = FALSE),
        error = function(e) NA_character_)
    if (length(out) == 0) "false" else "true"
}

#' Deterministic seed derived from run identity (AGENTS.md 9).
#'
#' "Use deterministic seeds derived from run ID, region, VMR task, repeat, and
#' fold." Same inputs always give the same seed; different inputs essentially
#' never collide.
seed_for <- function(run_id, region = "", task = "", repeat_i = "", fold = "") {
    key <- paste(run_id, region, task, repeat_i, fold, sep = "|")
    if (requireNamespace("digest", quietly = TRUE)) {
        h <- digest::digest(key, algo = "xxhash32", serialize = FALSE)
        return(strtoi(substr(h, 1, 7), base = 16L))
    }
    ## Fallback: sum of character codes with position weights. Weaker, but
    ## deterministic, which is the property that matters.
    codes <- utf8ToInt(key)
    as.integer(sum(codes * seq_along(codes)) %% .Machine$integer.max)
}

#' Write a table atomically: temp file in the same directory, then rename.
#'
#' A killed job leaves either the old file or no file, never a half-written one
#' that a downstream step would happily read.
write_atomic <- function(x, path, sep = "\t", ...) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(path, ".tmp.", Sys.getpid())
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    if (is.data.frame(x) || data.table::is.data.table(x)) {
        data.table::fwrite(x, tmp, sep = sep, ...)
    } else {
        writeLines(as.character(x), tmp)
    }
    if (!file.rename(tmp, path)) {
        stop("Atomic rename failed: ", tmp, " -> ", path)
    }
    invisible(path)
}

#' Reconcile a fanned-out set of tasks and refuse to close on unexplained loss.
#'
#' Every expected task must land in exactly one of: completed, excluded (by
#' policy, e.g. sex chromosomes), qc_failed (by a documented QC rule), or
#' failed (computational). Anything left over is unaccounted for.
#'
#' AGENTS.md 9: production runs have zero tolerance for unexplained
#' computational failures, so a nonzero `failed` count stops the run unless the
#' caller explicitly allows it (smoke runs).
reconcile <- function(expected, completed, excluded = character(),
                      qc_failed = character(), failed = character(),
                      run = NULL, allow_failures = FALSE) {
    expected  <- as.character(expected)
    completed <- as.character(completed)
    excluded  <- as.character(excluded)
    qc_failed <- as.character(qc_failed)
    failed    <- as.character(failed)

    accounted <- c(completed, excluded, qc_failed, failed)
    assert_no_dups(accounted, "reconciled task IDs (a task is in two categories)")

    unaccounted <- setdiff(expected, accounted)
    unexpected  <- setdiff(accounted, expected)

    summary <- data.table::data.table(
        category = c("expected", "completed", "excluded", "qc_failed",
                     "failed", "unaccounted", "unexpected"),
        n = c(length(expected), length(completed), length(excluded),
              length(qc_failed), length(failed), length(unaccounted),
              length(unexpected))
    )
    print(summary)

    if (!is.null(run)) {
        write_atomic(summary, file.path(run$dir, "task_reconciliation.tsv"))
        if (length(unaccounted) > 0 || length(failed) > 0) {
            write_atomic(
                data.table::data.table(
                    task = c(unaccounted, failed),
                    status = c(rep("unaccounted", length(unaccounted)),
                               rep("failed", length(failed)))),
                file.path(run$dir, "task_failures.tsv"))
        }
    }

    if (length(unexpected) > 0) {
        stop("Reconciliation found ", length(unexpected), " task(s) not in the ",
             "expected set: ", paste(head(unexpected, 10), collapse = ", "))
    }
    if (length(unaccounted) > 0) {
        stop("Reconciliation found ", length(unaccounted), " unaccounted task(s): ",
             paste(head(unaccounted, 10), collapse = ", "),
             "\n  Every expected task must be completed, excluded, QC-failed, ",
             "or failed. A missing output file is not an explanation.")
    }
    if (length(failed) > 0 && !allow_failures) {
        stop(length(failed), " task(s) failed computationally: ",
             paste(head(failed, 10), collapse = ", "),
             "\n  AGENTS.md 9: production runs have zero tolerance for ",
             "unexplained computational failures. Fix and rerun.")
    }
    invisible(summary)
}

#' Checksum every file in a run directory and close the run.
#'
#' After this, the directory is finished. Downstream modules cite the run ID.
close_run <- function(run, outputs = NULL) {
    if (is.null(outputs)) {
        outputs <- list.files(run$dir, recursive = TRUE, full.names = TRUE)
        outputs <- outputs[!grepl("(manifest\\.tsv|output_checksums\\.tsv)$", outputs)]
    }
    sums <- data.table::data.table(
        file = sub(paste0("^", run$dir, "/"), "", outputs),
        sha256 = vapply(outputs, file_sha256, character(1)),
        bytes = vapply(outputs, function(f) as.character(file.info(f)$size), character(1))
    )
    write_atomic(sums, file.path(run$dir, "output_checksums.tsv"))
    append_manifest(run, list(
        finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        n_output_files = nrow(sums)
    ))
    ## Make the run read-only. Immutability enforced by the filesystem, not by
    ## everyone remembering the rule.
    Sys.chmod(list.files(run$dir, recursive = TRUE, full.names = TRUE), mode = "0444")
    message("[run] closed ", run$dir, " (", nrow(sums), " output files)")
    invisible(sums)
}
