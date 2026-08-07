suppressPackageStartupMessages({
  library(data.table)
  library(quadprog)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

normalize_region <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z[z %in% c("dorsolateral prefrontal cortex", "prefrontal cortex")] <- "dlpfc"
  z[z %in% c("caudate nucleus")] <- "caudate"
  z
}

load_reference_objects <- function(path, prefix, platform = c("WGBS", "450850")) {
  platform <- match.arg(platform)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  suffix <- if (platform == "WGBS") "WGBS" else "450850"
  signature_name <- if (platform == "WGBS") {
    paste0(prefix, "_sig_all_WGBS")
  } else {
    paste0(prefix, "_sig_all")
  }
  list(
    dm = env[[paste0(prefix, "_DF_", suffix)]],
    signature = env[[signature_name]]
  )
}

top_reference_markers <- function(reference, prefix, n_per_cell) {
  dm <- reference$dm
  sig <- reference$signature
  dm <- dm[dm$TargetID %in% rownames(sig), , drop = FALSE]
  rbindlist(lapply(colnames(sig), function(cell_type) {
    ord <- order(dm[[cell_type]], na.last = NA)
    take <- head(ord, n_per_cell)
    data.table(
      reference = prefix,
      cell_type = cell_type,
      target_id = as.character(dm$TargetID[take]),
      marker_rank = seq_along(take),
      marker_pvalue = as.numeric(dm[[cell_type]][take])
    )
  }))
}

parse_reference_coordinates <- function(target_id) {
  clean <- sub("^chr", "", as.character(target_id))
  fields <- tstrsplit(clean, ":", fixed = TRUE)
  data.table(
    seqnames_hg19 = paste0("chr", fields[[1]]),
    pos_hg19 = suppressWarnings(as.integer(fields[[2]]))
  )
}

normalize_fraction_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  totals <- rowSums(x)
  zero <- !is.finite(totals) | totals <= 0
  if (any(zero)) x[zero, ] <- 1 / ncol(x)
  x / rowSums(x)
}

orient_fraction_matrix <- function(x, sample_ids, cell_types) {
  x <- as.matrix(x)
  if (all(sample_ids %in% rownames(x))) {
    x <- x[sample_ids, , drop = FALSE]
  } else if (all(sample_ids %in% colnames(x))) {
    x <- t(x[, sample_ids, drop = FALSE])
  } else if (nrow(x) == length(sample_ids)) {
    rownames(x) <- sample_ids
  } else if (ncol(x) == length(sample_ids)) {
    x <- t(x)
    rownames(x) <- sample_ids
  } else {
    stop("Cannot orient fraction matrix to bulk samples")
  }
  if (!all(cell_types %in% colnames(x))) {
    stop("Fraction matrix is missing expected cell types: ",
         paste(setdiff(cell_types, colnames(x)), collapse = ", "))
  }
  normalize_fraction_matrix(x[, cell_types, drop = FALSE])
}

build_reference_signature <- function(reference, bulk_ids, n_per_cell) {
  dm <- reference$dm
  sig <- reference$signature
  available <- intersect(intersect(rownames(sig), as.character(dm$TargetID)), bulk_ids)
  dm <- dm[dm$TargetID %in% available, , drop = FALSE]
  markers <- unique(unlist(lapply(colnames(sig), function(cell_type) {
    head(as.character(dm$TargetID[order(dm[[cell_type]], na.last = NA)]), n_per_cell)
  })))
  sig[intersect(markers, bulk_ids), , drop = FALSE]
}

simplex_qp_deconvolution <- function(bulk, signature, ridge = 1e-8) {
  common <- intersect(rownames(signature), rownames(bulk))
  ref <- as.matrix(signature[common, , drop = FALSE])
  obs <- as.matrix(bulk[common, , drop = FALSE])
  k <- ncol(ref)
  dmat <- 2 * crossprod(ref) + diag(ridge, k)
  amat <- cbind(rep(1, k), diag(k))
  bvec <- c(1, rep(0, k))
  out <- t(vapply(seq_len(ncol(obs)), function(j) {
    fit <- solve.QP(
      Dmat = dmat,
      dvec = 2 * crossprod(ref, obs[, j]),
      Amat = amat,
      bvec = bvec,
      meq = 1
    )
    pmax(0, fit$solution)
  }, numeric(k)))
  rownames(out) <- colnames(obs)
  colnames(out) <- colnames(ref)
  normalize_fraction_matrix(out)
}

run_reference_ensemble <- function(bulk, references, n_per_cell) {
  if (!requireNamespace("EpiDISH", quietly = TRUE)) {
    stop("EpiDISH is required for the recorded reference-ensemble fallback")
  }
  sample_ids <- colnames(bulk)
  components <- list()
  signature_summary <- list()
  for (prefix in names(references)) {
    sig <- build_reference_signature(references[[prefix]], rownames(bulk), n_per_cell)
    b <- bulk[rownames(sig), , drop = FALSE]
    cell_types <- colnames(sig)
    signature_summary[[prefix]] <- data.table(
      reference = prefix,
      n_signature_markers = nrow(sig),
      n_cell_types = ncol(sig)
    )
    for (method in c("RPC", "CP")) {
      invisible(capture.output(
        fit <- suppressMessages(
          tryCatch(
            EpiDISH::epidish(beta.m = b, ref.m = sig, method = method)$estF,
            error = function(e) e
          )
        )
      ))
      if (!inherits(fit, "error")) {
        components[[paste(prefix, method, sep = "_")]] <-
          orient_fraction_matrix(fit, sample_ids, cell_types)
      } else {
        warning(prefix, " ", method, " failed: ", conditionMessage(fit))
      }
    }
    components[[paste(prefix, "simplex_QP", sep = "_")]] <-
      simplex_qp_deconvolution(b, sig)
  }
  if (length(components) < 2) stop("Fewer than two deconvolution components succeeded")
  common_cells <- Reduce(intersect, lapply(components, colnames))
  components <- lapply(components, function(x) x[sample_ids, common_cells, drop = FALSE])
  arr <- simplify2array(components)
  ensemble <- apply(arr, c(1, 2), median, na.rm = TRUE)
  dimnames(ensemble) <- list(sample_ids, common_cells)
  list(
    ensemble = normalize_fraction_matrix(ensemble),
    components = components,
    signature_summary = rbindlist(signature_summary)
  )
}

standardize_cell_type_names <- function(x, label_map) {
  out <- unname(vapply(as.character(x), function(z) label_map[[z]] %||% z, character(1)))
  out
}
