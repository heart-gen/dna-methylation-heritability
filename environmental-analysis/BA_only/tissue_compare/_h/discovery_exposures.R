DISCOVERY_BINARY_ENVS <- c(
  "smoking",
  "codeine",
  "morphine",
  "cocaine",
  "ethanol",
  "antipsychotics",
  "nicotine",
  "amphetamines",
  "any_trauma_hx"
)

DISCOVERY_CATEGORICAL_ENVS <- c(
  "education",
  "marital_status"
)

get_discovery_env_vars <- function() {
  c(DISCOVERY_BINARY_ENVS, DISCOVERY_CATEGORICAL_ENVS)
}

report_env_var_coverage <- function(
  pheno_matrix,
  env_vars = get_discovery_env_vars()
) {
  sample_df <- pheno_matrix %>%
    dplyr::distinct(brnum, .keep_all = TRUE)

  coverage_rows <- lapply(env_vars, function(env) {
    valid <- sample_df %>%
      dplyr::filter(!is.na(.data[[env]]))

    counts <- table(valid[[env]], useNA = "no")
    count_str <- if (length(counts) == 0) {
      "none"
    } else {
      paste(
        sprintf("%s=%s", names(counts), as.integer(counts)),
        collapse = ", "
      )
    }

    data.frame(
      variable = env,
      n_total = nrow(sample_df),
      n_non_missing = nrow(valid),
      class_counts = count_str,
      stringsAsFactors = FALSE
    )
  })

  coverage <- dplyr::bind_rows(coverage_rows)
  print(coverage)
  invisible(coverage)
}
