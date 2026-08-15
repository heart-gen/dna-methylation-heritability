#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    input = file.path(dirname(script_path), "..", "_m", "observed"),
    expected = "",
    output_dir = file.path(dirname(script_path), "..", "_m", "observed", "combined"),
    max_failure_rate = "0.05",
    min_within_domain_rate = "0.80",
    fail_on_qc = "TRUE"
))
if (!nzchar(cli$expected)) stop("--expected is required")
expected <- read_tsv(cli$expected)
required_expected <- c("region", "population", "expected_tasks")
if (!all(required_expected %in% names(expected))) {
    stop("Expected-task manifest must contain: ", paste(required_expected, collapse = ", "))
}
read_category <- function(category) {
    files <- list.files(
        cli$input,
        pattern = "vmr-[0-9]+\\.tsv$",
        full.names = TRUE,
        recursive = TRUE
    )
    files <- files[grepl(paste0("/", category, "/"), files)]
    if (!length(files)) return(NULL)
    tables <- lapply(files, read_tsv)
    columns <- Reduce(union, lapply(tables, names))
    tables <- lapply(tables, function(data) {
        missing <- setdiff(columns, names(data))
        for (column in missing) data[[column]] <- NA
        data[, columns, drop = FALSE]
    })
    do.call(rbind, tables)
}
summaries <- read_category("summary")
failures <- read_category("failures")
qc_failures <- read_category("qc_failures")
excluded <- read_category("excluded")
if (is.null(summaries)) stop("No observed VMR estimates were produced")

summaries <- summaries[order(
    summaries$region, summaries$population,
    as.character(summaries$chromosome), summaries$start
), , drop = FALSE]
if (anyDuplicated(paste(summaries$region, summaries$population, summaries$task_id))) {
    stop("Duplicate analyzed task IDs detected")
}
task_keys <- function(data) {
    if (is.null(data)) return(character())
    paste(data$region, data$population, data$task_id)
}
all_keys <- c(
    task_keys(summaries), task_keys(excluded), task_keys(qc_failures),
    task_keys(failures)
)
if (anyDuplicated(all_keys)) {
    stop("A task ID was recorded in more than one terminal category")
}

count_category <- function(data, region, population) {
    if (is.null(data)) return(0L)
    sum(data$region == region & data$population == population)
}
qc <- lapply(seq_len(nrow(expected)), function(i) {
    region <- expected$region[[i]]
    population <- expected$population[[i]]
    analyzed <- count_category(summaries, region, population)
    computational_failed <- count_category(failures, region, population)
    qc_failed <- count_category(qc_failures, region, population)
    excluded_count <- count_category(excluded, region, population)
    region_summary <- summaries[
        summaries$region == region & summaries$population == population,
        , drop = FALSE
    ]
    failed <- computational_failed + qc_failed
    recorded <- analyzed + failed + excluded_count
    failure_rate <- failed / expected$expected_tasks[[i]]
    qc_failure_rate <- qc_failed / expected$expected_tasks[[i]]
    computational_failure_rate <- computational_failed / expected$expected_tasks[[i]]
    within_domain_rate <- if (analyzed) {
        mean(region_summary$calibration_status == "within_domain")
    } else {
        NA_real_
    }
    data.frame(
        region = region,
        population = population,
        expected_tasks = expected$expected_tasks[[i]],
        recorded_tasks = recorded,
        analyzed_tasks = analyzed,
        excluded_tasks = excluded_count,
        failed_tasks = failed,
        failure_rate = failure_rate,
        qc_failed_tasks = qc_failed,
        qc_failure_rate = qc_failure_rate,
        computational_failed_tasks = computational_failed,
        computational_failure_rate = computational_failure_rate,
        within_domain_rate = within_domain_rate,
        positive_signal_rate = if (analyzed) mean(region_summary$positive_signal) else NA_real_,
        complete = recorded == expected$expected_tasks[[i]],
        failure_qc_pass = qc_failure_rate <=
            as_num(cli$max_failure_rate, "max_failure_rate"),
        computational_qc_pass = computational_failed == 0L,
        domain_qc_pass = is.finite(within_domain_rate) &&
            within_domain_rate >= as_num(cli$min_within_domain_rate, "min_within_domain_rate"),
        stringsAsFactors = FALSE
    )
})
qc <- do.call(rbind, qc)
qc$overall_qc_pass <- qc$complete & qc$failure_qc_pass &
    qc$computational_qc_pass & qc$domain_qc_pass

dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(summaries, file.path(cli$output_dir, "calibrated-local-h2-AA-vmrs.tsv"))
if (!is.null(failures)) write_tsv(failures, file.path(cli$output_dir, "failed-vmrs.tsv"))
if (!is.null(qc_failures)) {
    write_tsv(qc_failures, file.path(cli$output_dir, "qc-failed-vmrs.tsv"))
}
if (!is.null(excluded)) write_tsv(excluded, file.path(cli$output_dir, "excluded-vmrs.tsv"))
write_tsv(qc, file.path(cli$output_dir, "observed-run-qc.tsv"))
status <- as.data.frame(table(
    region = summaries$region,
    population = summaries$population,
    calibration_status = summaries$calibration_status
), stringsAsFactors = FALSE)
write_tsv(status, file.path(cli$output_dir, "calibration-status-counts.tsv"))
capture_session_info(file.path(cli$output_dir, "session-info.txt"))
cat("Combined", nrow(summaries), "observed VMR estimates\n")
if (as_bool(cli$fail_on_qc, "fail_on_qc") && !all(qc$overall_qc_pass)) {
    quit(save = "no", status = 2L)
}
