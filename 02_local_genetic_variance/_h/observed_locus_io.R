## Shim. The observed-locus reader was promoted to 00_shared/locus_io.R on
## 2026-08-23 so that Module 03 consumes byte-identical genotype, covariate and
## donor-alignment logic -- 02 and 03 must score the same variants and donors or
## their endpoints are not comparable.
##
## The function body did not change. This file exists so Module 02's Stage 01,
## the observed-regime scenario runner, and the geometry scanner keep working
## unedited; those scripts are part of six accepted, sealed production runs.
source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "locus_io.R"))
