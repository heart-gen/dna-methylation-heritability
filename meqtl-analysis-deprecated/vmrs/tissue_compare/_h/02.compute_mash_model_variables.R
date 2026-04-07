#############################################
## Mash modeling in R for xQTL
#############################################
suppressPackageStartupMessages({
    library(mashr)
    library(dplyr)
    library(ggpubr)
})

save_img <- function(image, output_dir, fn, w = 7, h = 7) {
    for (ext in c(".svg", ".pdf")) {
        ggsave(file = file.path(output_dir, paste0(fn, ext)),
               plot = image, width = w, height = h)
    }
}

get_bhat <- function(input_dir, file = "bhat_nominal_3regions_AA.txt.gz") {
    data.table::fread(file.path(input_dir, file), header = TRUE, sep = "\t") |>
        mutate(effect = paste(shared_feature_id, variant_id, sep = "_")) |>
        distinct(effect, .keep_all = TRUE)
}

get_shat <- function(input_dir, file = "shat_nominal_3regions_AA.txt.gz") {
    data.table::fread(file.path(input_dir, file), header = TRUE, sep = "\t") |>
        mutate(effect = paste(shared_feature_id, variant_id, sep = "_")) |>
        distinct(effect, .keep_all = TRUE)
}

plot_mixture_prop <- function(m, output_dir) {
    fn <- "barplot_estimated_pi"
    df <- get_estimated_pi(m) |> as.data.frame() |>
        tibble::rownames_to_column("Model")

    colnames(df)[2] <- "Estimated pi"
    brp <- ggbarplot(
        df, x = "Model", y = "Estimated pi", fill = "gray",
        ggtheme = theme_pubr(base_size = 20), xlab = "",
        label = TRUE, label.pos = "out", lab.nb.digits = 2
    ) + font("y.title", face = "bold") + rotate_x_text(45)

    save_img(brp, output_dir, fn, 7, 7)
}

run_mashr <- function(seed, percentage, input_dir, output_dir) {
    set.seed(seed)
                                        # Load prepared dat
    message("Load prepared data")
    bhat <- get_bhat(input_dir) |>
        tibble::column_to_rownames("effect") |>
        select(-shared_feature_id, -variant_id) |>
        as.matrix()

    shat <- get_shat(input_dir) |>
        tibble::column_to_rownames("effect") |>
        select(-shared_feature_id, -variant_id) |>
        as.matrix()

    save(bhat, shat, file = file.path(output_dir, "bhat_shat.RData"))
    n_conditions <- ncol(bhat)

                                        # Prepare random and strong subsets
    message("Prepare random and strong subsets")
    rand_n <- max(2L, round(nrow(bhat) * percentage))

    if (rand_n >= nrow(bhat)) {
        random_subset <- seq_len(nrow(bhat))
    } else {
        random_subset <- sample.int(nrow(bhat), rand_n)
    }

    m_1by1 <- mash_1by1(mash_set_data(bhat, shat))
    strong_subset <- get_significant_results(m_1by1, 0.05)

    if (!length(strong_subset)) {
        stop("No strong effects detected; cannot fit data-driven covariances.")
    }

    rm(m_1by1)
    gc(verbose = TRUE)

                                        # Correlation structure
    message("Estimate correlation structure")
    data_init <- mash_set_data(bhat[random_subset, ], shat[random_subset, ])
    Vhat <- estimate_null_correlation_simple(data_init)
    rm(data_init)
    gc(verbose = TRUE)

                                        # Apply correlation structure
    message("Apply correlation structure to random and strong subsets")
    data_random <- mash_set_data(bhat[random_subset, ],
                                 shat[random_subset, ], V = Vhat)
    data_strong <- mash_set_data(bhat[strong_subset, ],
                                 shat[strong_subset, ], V = Vhat)
    rm(bhat, shat)
    gc(verbose = TRUE)

                                        # Data-driven covariances
    message("Compute data-driven covariances")
    U_pca <- cov_pca(data_strong, npc = n_conditions)
    U_ed  <- cov_ed(data_strong, U_pca)
    rm(U_pca)
    gc(verbose = TRUE)
                                        # Fit mash model
    message("Fit mash model")
    U_c <- cov_canonical(data_random)
    m   <- mash(data_random, Ulist = c(U_ed, U_c), outputlevel = 1)
    rm(U_ed, U_c, data_random)
    gc(verbose = TRUE)

                                        # Posterior summaries
    message("Compute posterior summaries")
    m2 <- mash(data_strong, g = get_fitted_g(m), fixg = TRUE)
    save(m, Vhat, file = file.path(output_dir, "model_variables.RData"))
    rm(data_strong, m, Vhat)
    gc(verbose = TRUE)

                                        # Pairwise sharing
    message("Pairwise sharing")
    print(get_pairwise_sharing(m2))

                                        # Significant results
    message("Save significant results")
    sig_idx <- get_significant_results(m2)
    print(length(sig_idx))
    save(m2, file = file.path(output_dir, "mashr_meta_results.RData"))

    # Export LFSR and posterior means
    message("Export LFSR + posterior means")
    m2$result$lfsr |> as.data.frame() |>
        tibble::rownames_to_column("effect") |>
        data.table::fwrite(file.path(output_dir, "lfsr.strong_signals.tsv"),
                           sep = "\t")
    m2$result$PosteriorMean |> as.data.frame() |>
        tibble::rownames_to_column("effect") |>
        data.table::fwrite(file.path(output_dir, "posterior_mean.strong_signals.tsv"),
                           sep = "\t")

                                        # Mixture proportions
    message("Estimate mixture proportions")
    print(get_estimated_pi(m2))
    plot_mixture_prop(m2, output_dir)
}

## Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("
Usage: Rscript 02.compute_mash_model_variables.R <input_dir> [output_dir] [seed] [percentage]

Required:
  input_dir       Directory containing bhat/shat files

Optional:
  output_dir      Where results go (default = input_dir)
  seed            Default = 20220422
  percentage      Random subset percentage, default = 0.05
")
}

input_dir   <- normalizePath(args[1], mustWork = TRUE)
output_dir  <- ifelse(length(args) >= 2, normalizePath(args[2], mustWork = FALSE), input_dir)
seed        <- ifelse(length(args) >= 3, as.integer(args[3]), 20220422L)
percentage  <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.05)

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

run_mashr(seed, percentage, input_dir, output_dir)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
