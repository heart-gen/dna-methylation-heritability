## Run a specific chunk with mashr (no argparse)
suppressPackageStartupMessages({
  library(mashr)
  library(dplyr)
})

usage <- function() {
    cat(
        "Usage: Rscript 04.allpairs_mash_model.R -c <chunk_num> [-i <model_dir>] [-o <output_dir>]\n",
        "  -c, --chunk       Integer chunk number to run (required)\n",
        "  -m, --model       Model directory with RData files [default: model]\n",
        "  -o, --output      Output directory for results    [default: output]\n",
        "  -h, --help        Show this help and exit\n",
        sep = ""
    )
}

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    out <- list(chunk = NA_integer_, model = "model", output = "output")

    i <- 1L
    next_val <- function() { # Parse flag
        i <<- i + 1L
        if (i > length(args)) stop("Missing value after flag.", call. = FALSE)
        args[[i]]
    }

    while (i <= length(args)) {
        a <- args[[i]]
        if (a %in% c("-h", "--help")) {
            usage()
            quit(save = "no", status = 0)
        } else if (grepl("^--chunk=", a)) {
            out$chunk <- as.integer(sub("^--chunk=", "", a))
        } else if (a %in% c("-c", "--chunk")) {
            out$chunk <- as.integer(next_val())
        } else if (grepl("^--model=", a)) {
            out$model <- sub("^--model=", "", a)
        } else if (a %in% c("-m", "--model")) {
            out$model <- next_val()
        } else if (grepl("^--output=", a)) {
            out$output <- sub("^--output=", "", a)
        } else if (a %in% c("-o", "--output")) {
            out$output <- next_val()
        } else {
            stop(sprintf("Unknown argument: %s", a), call. = FALSE)
        }
        i <- i + 1L
    }

    if (is.na(out$chunk)) {
        usage()
        stop("`--chunk` is required.", call. = FALSE)
    }

    out
}

ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_results_chunk <- function(chunk_num, model_dir, outdir) {
                                        # Load model + chunk data
    load(file.path(model_dir, "model_variables.RData"))
    load(file.path(outdir, paste0("chunk_", chunk_num, "_bhat_shat.RData")))

                                        # Run mash on the chunk
    data_chunk <- mash_set_data(bhat_chunk, shat_chunk, V = Vhat)
    m_chunk    <- mash(data_chunk, g = get_fitted_g(m), fixg = TRUE)

                                        # lfsr
    fn_lfsr <- file.path(outdir, paste0("lfsr_", chunk_num, "_ancestry.txt.gz"))
    lfsr_df <- data.frame(m_chunk$result$lfsr) |>
        tibble::rownames_to_column("effect") |>
        tidyr::separate(col = effect, into = c("phenotype_id", "variant_id"),
                        remove = FALSE, sep = "_", extra = "merge"
                        )
    readr::write_tsv(lfsr_df, fn_lfsr)

                                        # posterior means (drop suffix `.mean` if present)
    fn_post <- file.path(outdir, paste0("posterior_mean_", chunk_num, "_ancestry.txt.gz"))
    post_df <- data.frame(m_chunk$result$PosteriorMean) |>
        dplyr::rename_with(~ stringr::str_replace_all(.x, "\\.mean$", ""),
                           .cols = dplyr::everything()) |>
        tibble::rownames_to_column("effect") |>
        tidyr::separate(col = effect, into = c("phenotype_id", "variant_id"),
                        remove = FALSE, sep = "_", extra = "merge"
                        )
    readr::write_tsv(post_df, fn_post)
    invisible(TRUE)
}

## MAIN
args <- parse_args()
ensure_dir(args$output)

message(sprintf("Run chunk: %s", args$chunk))
message(sprintf("Model dir: %s", normalizePath(args$model, mustWork = FALSE)))
message(sprintf("Output dir: %s", normalizePath(args$output, mustWork = FALSE)))

save_results_chunk(args$chunk, args$model, args$output)

# Reproducibility info
print(Sys.time())
print(proc.time())
options(width = 120)
sessioninfo::session_info()
