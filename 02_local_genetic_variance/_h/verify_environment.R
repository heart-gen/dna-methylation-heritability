#!/usr/bin/env Rscript

required <- c(
    glmnet = "4.1.0",
    bigsnpr = "1.12.0",
    ggplot2 = "3.5.0",
    sessioninfo = "1.2.0"
)
missing <- names(required)[!vapply(
    names(required), requireNamespace, logical(1L), quietly = TRUE
)]
if (length(missing)) {
    stop("Required R packages are unavailable: ", paste(missing, collapse = ", "))
}
too_old <- names(required)[vapply(names(required), function(package) {
    utils::packageVersion(package) < numeric_version(required[[package]])
}, logical(1L))]
if (length(too_old)) {
    stop("R packages below required versions: ", paste(too_old, collapse = ", "))
}
if (getRversion() < "4.3.0" || getRversion() >= "4.5.0") {
    stop("R must be >=4.3.0 and <4.5.0; found ", getRversion())
}
cat("Environment verification passed\n")
cat("R:", as.character(getRversion()), "\n")
for (package in names(required)) {
    cat(package, ":", as.character(utils::packageVersion(package)), "\n")
}
