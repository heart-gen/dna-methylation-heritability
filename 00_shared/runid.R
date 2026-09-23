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

    commit <- git_commit(root)
    has_code <- git_commit_has_module_code(commit, module_root, root)
    if (identical(has_code, "false")) {
        ## Loud, because a SLURM log is where this will be read. NOT fatal:
        ## refusing the run is a policy change (AGENTS.md 12), and the case this
        ## catches is a module whose code is not yet committed at all -- exactly
        ## 01b_estimation_cells' first run. The field below makes the defect
        ## machine-detectable in the sealed manifest either way; promoting it to
        ## a stop() is a one-line change for the PI to authorise.
        msg <- paste0(
            "Recorded git_commit ", commit, " contains no files under ",
            module_root, ".\n",
            "  This run will not be reproducible from the commit it records ",
            "(AGENTS.md 9).\n",
            "  Commit the module's _h/ code before a production run.")
        warning(msg, call. = FALSE, immediate. = TRUE)
        message("[run] PROVENANCE WARNING: ", msg)
    }

    manifest <- c(
        list(
            run_id          = run_id,
            analysis        = module,
            cohort          = cohort,
            region          = region,
            vmr_set_id      = vmr_set_id,
            git_commit      = commit,
            git_dirty       = git_dirty(root),
            git_commit_has_module_code = has_code,
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

#' Is a run sealed? TRUE only when the manifest records a non-empty finished_at.
#'
#' Absent row, NA and "" all mean the run is still open: seal_run() is the only
#' writer of finished_at, so a freshly opened run has no such row at all.
#'
#' Written as a function because the obvious inline form is wrong in a way that
#' is invisible on a sealed run and fatal on an open one:
#' `manifest$value[manifest$field == "finished_at"][1]` on a manifest with no
#' such row yields NA_character_ (subsetting past the end, not a zero-length
#' vector), `%||%` replaces NULL but not NA, and nzchar(NA_character_) is TRUE.
#' Every stage guarded that way refuses to run on the run it just opened.
run_is_sealed <- function(manifest) {
    v <- manifest$value[manifest$field == "finished_at"]
    length(v) >= 1L && !is.na(v[[1L]]) && nzchar(v[[1L]])
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

## 00_shared/slurm.sh runs `module purge`, which strips Quest's git module from
## PATH. Compute nodes have no /usr/bin/git fallback, so a bare "git" here is not
## resolvable in a batch job. slurm.sh exports V2_GIT_BIN (resolved BEFORE the
## purge); fall back to "git" for interactive use.
git_bin <- function() {
    b <- Sys.getenv("V2_GIT_BIN", unset = "")
    if (nzchar(b) && file.exists(b)) b else "git"
}

## Returns NA when git could not be run at all, so a missing git is recorded as
## unknown rather than mistaken for a value.
git_run <- function(root, cmd) {
    tryCatch(
        suppressWarnings(system2(git_bin(), c("-C", shQuote(root), cmd),
                                 stdout = TRUE, stderr = FALSE)),
        error = function(e) NULL)
}

git_commit <- function(root = repo_root()) {
    out <- git_run(root, c("rev-parse", "HEAD"))
    if (is.null(out) || length(out) == 0) NA_character_ else out[1]
}

git_dirty <- function(root = repo_root()) {
    ## A FAILED git call previously returned character(0) and was reported as
    ## "false" -- a broken check claiming a clean tree. Distinguish the three
    ## states: unknown (git unusable), clean, dirty. `git status` exits 0 with
    ## empty output on a clean tree, so success is confirmed separately.
    ok <- tryCatch(
        suppressWarnings(system2(git_bin(), c("-C", shQuote(root), "rev-parse",
                                              "--is-inside-work-tree"),
                                 stdout = TRUE, stderr = FALSE)),
        error = function(e) NULL)
    if (is.null(ok) || length(ok) == 0 || !identical(ok[1], "true"))
        return(NA_character_)
    out <- git_run(root, c("status", "--porcelain"))
    if (is.null(out)) NA_character_ else if (length(out) == 0) "false" else "true"
}

#' Does `commit` actually contain the module's own code?
#'
#' AGENTS.md 9 requires every production output to carry a Git commit. A commit
#' that predates the module is worse than no commit: it looks like provenance
#' and reproduces nothing. Six sealed 01b_estimation_cells runs record
#' f8fd01a74, which contains zero 01b_estimation_cells/ files, so those runs
#' cannot be regenerated from what they record. Their configs and outputs are
#' intact -- this is a provenance defect, not a data defect -- but it is
#' mechanically checkable at run time, which is why it is checked here.
#'
#' Returns "true", "false", or NA_character_ when git could not answer. The
#' three states are kept distinct for the same reason git_dirty() keeps them:
#' a check that cannot run must not be recorded as a check that passed.
git_commit_has_module_code <- function(commit, module_root, root = repo_root()) {
    if (is.null(commit) || length(commit) == 0 || is.na(commit)) {
        return(NA_character_)
    }
    rel <- tryCatch({
        m <- normalizePath(module_root, mustWork = FALSE)
        r <- normalizePath(root, mustWork = FALSE)
        if (!startsWith(m, r)) return(NA_character_)
        sub("^/+", "", substring(m, nchar(r) + 1L))
    }, error = function(e) NA_character_)
    if (is.na(rel) || !nzchar(rel)) return(NA_character_)

    ## Confirm the commit resolves first. ls-tree on an unknown revision exits
    ## nonzero with empty stdout, which is indistinguishable from "the module is
    ## absent" unless the two questions are asked separately.
    ## shQuote because git_run() hands its args to system2(), which builds a
    ## /bin/sh command line rather than an argv: ^{commit} and the path would
    ## otherwise be exposed to the shell.
    ok <- git_run(root, c("rev-parse", "--verify", "--quiet",
                          shQuote(paste0(commit, "^{commit}"))))
    if (is.null(ok) || length(ok) == 0 || !nzchar(ok[1])) return(NA_character_)

    out <- git_run(root, c("ls-tree", "-r", "--name-only", shQuote(commit),
                           "--", shQuote(rel)))
    if (is.null(out)) return(NA_character_)
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
        ## Write the failure list UNCONDITIONALLY, header-only when there are no
        ## failures. It was previously written only when failures existed and was
        ## never cleared, so a run that failed tasks, was re-driven and then
        ## reconciled clean kept the old list beside a reconciliation reporting
        ## zero failures. A reader auditing AGENTS.md 9 compliance then sees
        ## failures that did not occur in the accepted attempt. Both files are
        ## now written by the same unconditional write_atomic() call in the same
        ## block, so the pair always describes one reconcile.
        write_atomic(
            data.table::data.table(
                task = c(unaccounted, failed),
                status = c(rep("unaccounted", length(unaccounted)),
                           rep("failed", length(failed)))),
            file.path(run$dir, "task_failures.tsv"))
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

#' Seal a run directory tree: clear every write bit, files AND directories.
#'
#' POSIX checks the *directory's* write bit when a name is created or removed,
#' never the file's own mode. A run whose files are all 0444 inside directories
#' that are still 0775 is therefore not immutable: any result in it can be
#' unlinked and written back under the same name, and new files can be added,
#' with no error raised and nothing to show for it in output_checksums.tsv.
#'
#' Until 2026-09-23 close_run() chmodded only files, so that was true of every
#' module except 02, which runs its own `chmod -R a-w` in
#' 06_finalize_observed_run.R afterwards. That asymmetry is how the defect was
#' found: during the run-retirement cleanup, Module 02 runs were the only ones
#' that resisted deletion, and every other module's "sealed" run deleted freely.
#' AGENTS.md 5.2 asks for immutability enforced by the filesystem rather than by
#' everyone remembering the rule, and half a seal does not do that.
#'
#' The mode is computed per directory as `mode & ~0222`, not set to a constant
#' "0555". These trees are setgid (2775) so that everything written beneath them
#' inherits the project group; a flat chmod to 0555 drops that bit silently, and
#' the loss only shows up later as a file owned by the wrong group. Clearing
#' just the write bits preserves setgid, and preserves the execute bits on files
#' that carry them, such as the `_h/` launchers Module 11 snapshots into `code/`.
#'
#' This cannot be applied retroactively. An already-sealed run's files are 0444
#' and its directories are whatever they were, and rewriting the mode of a run
#' that has been cited would change a run after it was closed, which is exactly
#' what 5.2 forbids. The fix binds new runs only; existing runs outside Module 02
#' keep writable directories and that is a fact about them, not a bug to repair.
#'
#' A seal that fails is worse than no seal, because everything downstream then
#' assumes an immutability that is not there -- this defect survived unnoticed
#' for exactly that reason. So verify and stop, rather than trusting the chmod.
seal_run_dir <- function(dir) {
    ## Without this, a typo'd path seals nothing and reports success: both
    ## listings below return character(0) on a directory that does not exist,
    ## and the verification pass then finds nothing writable.
    if (!dir.exists(dir)) stop("Cannot seal, no such directory: ", dir)

    ## all.files = TRUE because the default omits dotfiles, and a file this
    ## listing misses stays writable inside a run that claims to be immutable.
    ## That was the first of these seal defects (Module 11's _h/.gitkeep,
    ## 2026-09-22); with recursive = TRUE it adds hidden files only, never
    ## "." or "..".
    files <- list.files(dir, recursive = TRUE, full.names = TRUE,
                        all.files = TRUE)
    if (length(files) > 0) Sys.chmod(files, mode = "0444")

    ## list.dirs() returns `dir` itself and every subdirectory including hidden
    ## ones, which is what this needs: the top of the run must be sealed too,
    ## or its own entries stay unlinkable-and-replaceable.
    dirs <- list.dirs(dir, recursive = TRUE, full.names = TRUE)
    for (d in dirs) {
        mode <- file.info(d)$mode
        if (is.na(mode)) stop("Cannot stat directory while sealing: ", d)
        Sys.chmod(d, mode = as.octmode(bitwAnd(as.integer(mode),
                                               bitwNot(strtoi("222", 8L)))))
    }

    writable <- character(0)
    for (p in c(files, dirs)) {
        mode <- file.info(p)$mode
        if (!is.na(mode) && bitwAnd(as.integer(mode), strtoi("222", 8L)) != 0) {
            writable <- c(writable, p)
        }
    }
    if (length(writable) > 0) {
        stop("Seal incomplete: ", length(writable), " path(s) still writable ",
             "under ", dir, ", first few: ",
             paste(utils::head(writable, 5), collapse = ", "),
             "\n  The run is not immutable and must not be cited (AGENTS.md 5.2).")
    }
    invisible(list(files = length(files), dirs = length(dirs)))
}

#' Checksum every file in a run directory and close the run.
#'
#' After this, the directory is finished. Downstream modules cite the run ID.
close_run <- function(run, outputs = NULL) {
    ## Refuse to re-close a run that is already closed. Before the directory
    ## seal covered directories, this case SUCCEEDED: write_atomic() renames a
    ## temp file into place, rename replaces a 0444 target happily as long as
    ## the containing directory is writable, so re-running a finalizer quietly
    ## rewrote a sealed run's manifest and checksums -- an update in place, which
    ## AGENTS.md 5.2 forbids outright. Sealing the directories turns that into a
    ## permission error instead, which is better but reads like a filesystem
    ## problem. Say what it actually is, and what to do about it.
    mf <- file.path(run$dir, "manifest.tsv")
    if (file.exists(mf) &&
        run_is_sealed(data.table::fread(mf, colClasses = "character"))) {
        stop("Run is already closed: ", run$dir,
             "\n  Its manifest records a finished_at, so it has been sealed and ",
             "may already be cited.\n  AGENTS.md 5.2: never update a completed ",
             "run in place -- supersede it by minting a new run ID.")
    }

    if (is.null(outputs)) {
        ## all.files = TRUE, because the default omits dotfiles and a file this
        ## function does not list is neither checksummed nor made read-only --
        ## it stays writable inside a run that claims to be immutable. Found
        ## 2026-09-22, when Module 11 began snapshotting _h/ into the run and
        ## carried an _h/.gitkeep in with it. With recursive = TRUE this adds
        ## hidden files only: "." and ".." are not returned.
        outputs <- list.files(run$dir, recursive = TRUE, full.names = TRUE,
                              all.files = TRUE)
        ## The run's own two bookkeeping files, matched by exact path. This was
        ## a suffix pattern -- "(manifest\\.tsv|output_checksums\\.tsv)$" --
        ## which is unanchored at the front, so it also swallowed any OUTPUT
        ## whose name merely ends that way. Module 11's
        ## tables/software-and-run-manifest.tsv is one, and it went unchecksummed
        ## in fig-all-20260922-a without any warning. Name the two files; do not
        ## describe them.
        outputs <- setdiff(outputs,
                           file.path(run$dir, c("manifest.tsv",
                                                "output_checksums.tsv")))
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
    ## everyone remembering the rule. Directories are sealed as well as files;
    ## see seal_run_dir() for why chmodding only the files left every run
    ## outside Module 02 quietly mutable.
    sealed <- seal_run_dir(run$dir)
    message("[run] closed ", run$dir, " (", nrow(sums), " output files, ",
            "sealed ", sealed$files, " files and ", sealed$dirs, " directories)")
    invisible(sums)
}
