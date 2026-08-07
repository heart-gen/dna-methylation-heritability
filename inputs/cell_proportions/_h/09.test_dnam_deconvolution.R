#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
ROOT <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
source(file.path(ROOT, "inputs/cell_proportions/_h/dnam_deconvolution_utils.R"))

stopifnot(identical(normalize_region(c("Caudate nucleus", "DLPFC", "Hippocampus")),
                    c("caudate", "dlpfc", "hippocampus")))
coords <- parse_reference_coordinates(c("5:123", "chrX:456"))
stopifnot(identical(coords$seqnames_hg19, c("chr5", "chrX")))
stopifnot(identical(coords$pos_hg19, c(123L, 456L)))

set.seed(20260730)
signature <- matrix(runif(60), nrow = 20, ncol = 3,
                    dimnames = list(paste0("m", 1:20), c("Astro", "Micro", "Exc")))
truth <- matrix(c(0.2, 0.3, 0.5, 0.7, 0.1, 0.2), nrow = 2, byrow = TRUE,
                dimnames = list(c("S1", "S2"), colnames(signature)))
bulk <- signature %*% t(truth)
fit <- simplex_qp_deconvolution(bulk, signature)
stopifnot(max(abs(rowSums(fit) - 1)) < 1e-8)
stopifnot(max(abs(fit - truth)) < 1e-5)
stopifnot(all(fit >= 0 & fit <= 1))
cat("DNAm deconvolution utility tests passed\n")
