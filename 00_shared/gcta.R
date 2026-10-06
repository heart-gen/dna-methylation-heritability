#### GCTA / GREML I/O shared by 02b and 02c ####
##
## One implementation of the GCTA calls both cis-GREML modules make: the
## version pin, the REML call, the .hsq parser, the estimator-outcome
## classifier, and a GCTA-format GRM written straight from a dosage matrix.
## Moved here from 02b_greml_simulation_benchmark/_h/greml_functions.R on
## 2026-10-03 when 02c_cis_greml_sensitivity became the second caller
## (AGENTS.md 5.3: refactor shared logic into 00_shared, never duplicate it).
## Each function takes its settings as arguments; nothing here reads a config.

## ------------------------------------------------------------------- GCTA
gcta_version <- function(bin) {
    out <- suppressWarnings(system2(bin, character(0), stdout = TRUE,
                                    stderr = TRUE))
    hit <- regmatches(out, regexpr("v[0-9]+\\.[0-9]+(\\.[0-9]+)?", out))
    if (!length(hit)) stop("Could not read a version from ", bin)
    sub("^v", "", hit[[1]])
}

## Stop unless the binary is the pinned one. A REML result depends on the
## version (defaults and convergence criteria have changed across releases).
## `cfg` is any module config with a `gcta: {binary, version}` block.
assert_gcta_pin <- function(cfg) {
    bin <- cfg$gcta$binary
    if (!file.exists(bin)) stop("GCTA binary not found: ", bin)
    v <- gcta_version(bin)
    if (!identical(v, as.character(cfg$gcta$version))) {
        stop("GCTA version ", v, " does not match the configured pin ",
             cfg$gcta$version)
    }
    v
}

## Run gcta64 and classify the outcome.
##
## status is "completed" when an .hsq was written; "qc_failed" when GCTA stopped
## with one of the configured REML estimation messages (singular information
## matrix, non-convergence) -- an outcome of the estimator on that phenotype,
## counted and reported, not a crash; and "failed" for anything else, which
## reconcile() treats as an unexplained computational failure.
run_gcta <- function(bin, args, out_prefix, threads = 1L) {
    log <- paste0(out_prefix, ".gcta.log")
    st <- system2(bin, c(args, "--out", out_prefix, "--thread-num",
                         as.character(threads)),
                  stdout = log, stderr = log)
    invisible(list(exit = st, log = log))
}

classify_reml <- function(exit, hsq, log, failure_patterns, divergence = NULL) {
    if (identical(as.integer(exit), 0L) && file.exists(hsq)) {
        ## GCTA can write an .hsq whose h2 SE is 0 or absent: the information
        ## matrix was degenerate at the stopping point, so the estimate is not
        ## a REML optimum (02c saw values down to -128). Opt-in.
        if (isTRUE(divergence$degenerate_se_is_failure)) {
            se <- parse_hsq(hsq)$h2_se
            if (!is.finite(se) || se <= 0) {
                return(list(status = "qc_failed",
                            reason = "reml_estimation_failure: degenerate fit (h2 SE 0 or undefined)"))
            }
        }
        return(list(status = "completed", reason = NA_character_))
    }
    txt <- if (file.exists(log)) paste(readLines(log, warn = FALSE),
                                       collapse = "\n") else ""
    for (p in failure_patterns) {
        if (grepl(p, txt, fixed = TRUE)) {
            return(list(status = "qc_failed",
                        reason = paste0("reml_estimation_failure: ", p)))
        }
    }
    ## A diverging REML can crash GCTA instead of stopping it. Reconciled as an
    ## estimator outcome only when the exit code is a crash code AND GCTA's own
    ## iteration log shows the divergence; config/greml_benchmark.yml records
    ## the attempt that showed this.
    if (!is.null(divergence) &&
        as.integer(exit) %in% as.integer(unlist(divergence$exit_codes))) {
        nums <- suppressWarnings(as.numeric(
            regmatches(txt, gregexpr("-?[0-9]+\\.?[0-9]*(e[+-]?[0-9]+)?", txt))[[1]]))
        if (any(abs(nums) >= as.numeric(divergence$magnitude), na.rm = TRUE)) {
            return(list(status = "qc_failed",
                        reason = paste0("reml_estimation_failure: diverged, then gcta exit ",
                                        exit)))
        }
    }
    ## A diverging REML can also end in an ordinary GCTA error whose wording
    ## does not say so ("X^t * V^-1 * X is not invertible" also fires on a
    ## collinear covariate design). Such a message is an estimator outcome only
    ## when GCTA's own iteration trace shows a variance component that ran away
    ## from its EM starting value; a design fault stops before any iteration
    ## and so stays "failed". Opt-in: absent keys leave the rule off.
    gated <- unlist(divergence$gated_patterns)
    if (length(gated)) {
        hit <- gated[vapply(gated, grepl, logical(1), x = txt, fixed = TRUE)]
        ratio <- reml_trace_divergence(txt)
        if (length(hit) && is.finite(ratio) &&
            ratio >= as.numeric(divergence$relative_magnitude)) {
            return(list(status = "qc_failed",
                        reason = sprintf("reml_estimation_failure: diverged (max |V| %.2g x EM prior Vp), then: %s",
                                         ratio, hit[[1]])))
        }
    }
    tail_txt <- utils::tail(strsplit(txt, "\n")[[1]], 3)
    list(status = "failed",
         reason = paste0("gcta_exit_", exit, ": ",
                         paste(tail_txt, collapse = " | ")))
}

## Largest |variance component| in GCTA's REML iteration table, as a multiple
## of the summed EM-REML prior ("Updated prior values: ..."). NA when the log
## has no prior or no iteration rows, i.e. GCTA stopped before iterating.
reml_trace_divergence <- function(txt) {
    lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
    pl <- grep("^Updated prior values:", lines, value = TRUE)
    if (!length(pl)) return(NA_real_)
    prior <- sum(abs(suppressWarnings(as.numeric(
        strsplit(trimws(sub("^Updated prior values:", "", pl[[1]])), "[[:space:]]+")[[1]]))))
    hdr <- grep("^Iter\\.", lines)
    if (!length(hdr) || !is.finite(prior) || prior <= 0) return(NA_real_)
    it <- grep("^[0-9]+\t", lines[(hdr[[1]] + 1L):length(lines)], value = TRUE)
    if (!length(it)) return(NA_real_)
    v <- unlist(lapply(strsplit(it, "\t", fixed = TRUE), function(f)
        suppressWarnings(as.numeric(f[-(1:2)]))))
    if (!any(is.finite(v))) return(NA_real_)
    max(abs(v), na.rm = TRUE) / prior
}

## The estimator-outcome rules' settings, from config. The gated-pattern and
## degenerate-SE keys are optional; a config without them gets the crash-code
## rule only.
reml_divergence_rule <- function(cfg) {
    list(exit_codes = cfg$gcta$divergence_crash_exit_codes,
         magnitude = cfg$gcta$divergence_magnitude,
         gated_patterns = cfg$gcta$divergence_gated_patterns,
         relative_magnitude = cfg$gcta$divergence_relative_magnitude,
         degenerate_se_is_failure = cfg$gcta$degenerate_se_is_estimation_failure)
}

## Parse a GCTA .hsq into one row. Handles single-GRM ("V(G)/Vp") and
## multi-GRM ("Sum of V(G)/Vp") output; the multi-GRM total is what v1's
## 02.summary.py reported as `Sum of V(G)_Vp_Variance`.
parse_hsq <- function(path) {
    lines <- readLines(path, warn = FALSE)
    parts <- strsplit(lines[-1], "\t", fixed = TRUE)
    src <- vapply(parts, `[`, character(1), 1)
    val <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2)))
    se <- suppressWarnings(as.numeric(vapply(parts, function(x)
        if (length(x) >= 3) x[[3]] else NA_character_, character(1))))
    get <- function(s, v = val) { i <- match(s, src); if (is.na(i)) NA_real_ else v[[i]] }
    total <- if ("Sum of V(G)/Vp" %in% src) "Sum of V(G)/Vp" else "V(G)/Vp"
    data.table(
        h2_hat = get(total), h2_se = get(total, se),
        vg_hat = get("V(G)"), ve_hat = get("V(e)"), vp_hat = get("Vp"),
        logL = get("logL"), logL0 = get("logL0"), lrt = get("LRT"),
        pval = get("Pval"), n_used = get("n")
    )
}

## ------------------------------------------------------------------ the GRM
## A GCTA-format GRM written directly from a dosage matrix. Matches
## `gcta64 --make-grm` on the same SNPs and donors to float32 precision
## (max |diff| 1.0e-7 on a 118-donor, 3,428-SNP DLPFC locus with 1,930 missing
## calls; tests/test_greml_functions.R repeats that check). Missing calls
## contribute nothing and the denominator is the per-pair count of SNPs
## observed in both donors, which is GCTA's own rule. Writing it here rather
## than re-exporting a PLINK file keeps the genotype exactly the matrix that
## load_observed_locus() returned -- same window, same QC, same donor order.
grm_from_dosage <- function(G) {
    p <- colMeans(G, na.rm = TRUE) / 2
    keep <- is.finite(p) & p > 0 & p < 1
    G <- G[, keep, drop = FALSE]; p <- p[keep]
    Z <- sweep(G, 2, 2 * p, "-")
    Z <- sweep(Z, 2, sqrt(2 * p * (1 - p)), "/")
    M <- !is.na(Z)
    Z[!M] <- 0
    N <- tcrossprod(M * 1)
    list(A = tcrossprod(Z) / N, N = N, n_snps = ncol(G))
}

## GCTA stores the lower triangle row by row, which is the upper triangle in
## R's column-major order.
write_grm <- function(grm, ids, prefix) {
    A <- grm$A; N <- grm$N
    writeBin(as.numeric(A[upper.tri(A, diag = TRUE)]),
             paste0(prefix, ".grm.bin"), size = 4)
    writeBin(as.numeric(N[upper.tri(N, diag = TRUE)]),
             paste0(prefix, ".grm.N.bin"), size = 4)
    fwrite(as.data.table(ids)[, 1:2], paste0(prefix, ".grm.id"),
           sep = "\t", col.names = FALSE)
    invisible(prefix)
}

read_grm <- function(prefix) {
    ids <- fread(paste0(prefix, ".grm.id"), header = FALSE)
    n <- nrow(ids)
    v <- readBin(paste0(prefix, ".grm.bin"), "numeric", size = 4,
                 n = n * (n + 1) / 2)
    A <- matrix(0, n, n)
    A[upper.tri(A, diag = TRUE)] <- v
    A <- A + t(A) - diag(diag(A))
    list(A = A, ids = ids)
}

