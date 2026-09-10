#### 01b_estimation_cells / 00_new_run: open a donor-group estimation cell ####
##
## An ESTIMATION CELL is a donor subset of an accepted, sealed Module 01 run.
## The loci, the methylation phenotypes and the age/sex/diagnosis covariates all
## come from that pooled run unchanged; only the donors narrow. This stage picks
## the donors, proves the partition, and opens the immutable run directory that
## stages 01-04 fill.
##
## Why this is a separate module rather than a Module 01 option: Module 01 runs
## are immutable (AGENTS.md 5.2) and the pooled catalogs
## vmrcat-all_individuals-{region}-20260816 are accepted and SEALED. Per-group
## genotype extraction cannot be added to them in place, so it writes a new run
## that satisfies the same directory contract 00_shared/locus_io.R already
## consumes.
##
## This module NEVER re-derives VMRs, methylation residuals or phenotypes. It
## re-partitions donors over a fixed locus set. Anything else would make the two
## donor groups' loci non-identical, which is the entire point of the design
## (AGENTS.md 7.7: "Discovery happens once in the pooled sample").
##
## Usage:
##   Rscript 00_new_run.R --cohort all_individuals --group EA --region dlpfc \
##       --vmr-run-id vmrcat-all_individuals-dlpfc-20260816 [--allow-unlocked]

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

opts <- parse_v2_args(require = c("cohort", "region", "vmr_run_id"))

cell   <- opts$cell
region <- opts$region
group  <- opts$estimation_group

parsed <- parse_cell(cell)
if (!parsed$is_estimation_cell) {
    stop("01b materializes ESTIMATION CELLS. '", cell, "' is a discovery arm; ",
         "its Module 01 run is already the cell. Pass --group, e.g. ",
         "--cohort all_individuals --group EA.")
}

assert_locked(
    list(cohorts = load_config("cohorts"),
         thresholds = load_config("thresholds"),
         covariates_cells = list(
             pi_locked = load_config("covariates")$estimation_cells$pi_locked)),
    allow_unlocked = opts$allow_unlocked)

arm <- cell_def(cell)
catalog_cohort <- parsed$catalog_cohort

## ------------------------------------------------- upstream acceptance gate
##
## Module 01 predates the shared `Accepted runs` parser; its README is the
## record of acceptance. This is the same check
## 02_local_genetic_variance/_h/00_prepare_observed_run.R makes, kept identical
## on purpose -- 01b sits between 01 and 02 and must not be a weaker gate.
vmr_module_root <- file.path(V2_ROOT, "01_vmr_catalog")
vmr_run_dir <- file.path(vmr_module_root, "_m", "runs", opts$vmr_run_id)
if (!dir.exists(vmr_run_dir)) {
    stop("No such Module 01 run: ", vmr_run_dir)
}

vmr_readme <- readLines(file.path(vmr_module_root, "README.md"), warn = FALSE)
acceptance_row <- vmr_readme[grepl(paste0("| `", opts$vmr_run_id, "` |"),
                                   vmr_readme, fixed = TRUE)]
if (length(acceptance_row) != 1L ||
    !grepl("all five pass", acceptance_row, fixed = TRUE)) {
    stop("Module 01 run is not recorded as passing all five gates: ",
         opts$vmr_run_id)
}

vmr_manifest <- fread(file.path(vmr_run_dir, "manifest.tsv"),
                      colClasses = "character")
mval <- function(field) {
    v <- vmr_manifest$value[vmr_manifest$field == field]
    if (length(v) != 1L) stop("Module 01 manifest lacks unique field: ", field)
    as.character(v[[1L]])
}

## The catalog must be the cell's DISCOVERY arm, not the cell itself. This is
## the line that separates discovery from estimation: an EA cell is built on the
## POOLED catalog, never on an EA-only one, because an EA-only catalog would put
## the two donor groups on different loci.
if (!identical(mval("cohort"), catalog_cohort)) {
    stop("Cell '", cell, "' discovers on '", catalog_cohort,
         "' but Module 01 run ", opts$vmr_run_id, " is cohort '",
         mval("cohort"), "'")
}
if (!identical(tolower(mval("region")), region)) {
    stop("Requested region does not match Module 01 manifest")
}
if (!identical(toupper(mval("smoke_run")), "FALSE")) {
    stop("An estimation cell must start from a non-smoke Module 01 run")
}
vmr_set_id <- mval("vmr_set_id")

## --------------------------------------------------------- donor partition
##
## The pooled donor list is the ONLY source of candidate donors. Deriving the
## group from the phenotype table alone would let a donor who never entered the
## pooled catalog (no WGBS, no genotype) into the cell, so the cell would not be
## a subset of the discovery set and the two groups would not partition it.
pooled_donors <- fread(file.path(vmr_run_dir, "vmr", "donors_plink.txt"),
                       header = FALSE, colClasses = "character")
if (ncol(pooled_donors) < 2) {
    stop("Pooled donors_plink.txt has fewer than 2 columns")
}
setnames(pooled_donors, 1:2, c("FID", "IID"))
pooled_donors <- pooled_donors[, .(FID, IID)]
assert_no_dups(pooled_donors$FID, "FID in pooled donors_plink.txt")

pheno <- fread(arm$phenotype_table, header = TRUE)
for (col in c("brnum", "region", "race")) {
    if (!col %in% names(pheno)) {
        stop("Phenotype table is missing required column '", col, "': ",
             arm$phenotype_table)
    }
}
## Distinct name: `region` is also a column of the phenotype table, and the
## data.table `i` expression would resolve it to the column, not the argument.
## That is defect V2 in miniature.
region_arg <- region
pheno_region <- pheno[region == region_arg]
assert_no_dups(pheno_region$brnum, "brnum in phenotype table for this region")

## Every pooled donor must have a race label, or the two groups cannot be shown
## to partition the pooled set and a donor would vanish without a record.
race_of <- pheno_region$race[match(pooled_donors$FID, pheno_region$brnum)]
if (anyNA(race_of)) {
    missing <- pooled_donors$FID[is.na(race_of)]
    stop("Pooled donors absent from the phenotype table for ", region, ": ",
         paste(head(missing, 10), collapse = ", "))
}

cell_donors <- pooled_donors[race_of %in% arm$race_filter]
if (nrow(cell_donors) == 0) {
    stop("No pooled donors carry race in {", paste(arm$race_filter, collapse = ", "),
         "} for ", cell, "/", region)
}

## Partition proof. Every estimation cell on this catalog is checked, not just
## the one being built, so a mislabelled race_filter that made two cells overlap
## or lose a donor is caught here rather than in a downstream contrast.
cohorts_cfg <- load_config("cohorts")
sibling_cells <- Filter(
    function(nm) identical(cohorts_cfg$estimation_cells[[nm]]$catalog_cohort,
                           catalog_cohort),
    names(cohorts_cfg$estimation_cells))
sibling_filters <- lapply(sibling_cells, function(nm) {
    as.character(cohorts_cfg$estimation_cells[[nm]]$race_filter)
})
names(sibling_filters) <- sibling_cells

overlap <- unlist(sibling_filters)
if (anyDuplicated(overlap)) {
    stop("Estimation cells on '", catalog_cohort, "' share a race label: ",
         paste(unique(overlap[duplicated(overlap)]), collapse = ", "),
         ". The donor groups would not be disjoint.")
}
covered <- pooled_donors$FID[race_of %in% overlap]
uncovered <- setdiff(pooled_donors$FID, covered)
if (length(uncovered)) {
    stop("The estimation cells on '", catalog_cohort, "' do not cover every ",
         "pooled donor. Uncovered: ",
         paste(head(uncovered, 10), collapse = ", "),
         if (length(uncovered) > 10)
             paste0(" (and ", length(uncovered) - 10, " more)"))
}
if (nrow(cell_donors) >= nrow(pooled_donors)) {
    stop("Cell '", cell, "' is not a strict subset of the pooled donor set (",
         nrow(cell_donors), " of ", nrow(pooled_donors), ")")
}

assert_expected_n(nrow(cell_donors), cell, region)

## ------------------------------------------------------------ open the run
run <- new_run(
    module = "estcell",
    cohort = cell,
    region = region,
    module_root = file.path(V2_ROOT, "01b_estimation_cells"),
    vmr_set_id = vmr_set_id,
    upstream = list(vmr_catalog = opts$vmr_run_id),
    extra = list(
        smoke_run        = opts$allow_unlocked,
        catalog_cohort   = catalog_cohort,
        estimation_group = group,
        race_filter      = paste(arm$race_filter, collapse = ","),
        pgen_prefix      = arm$pgen_prefix,
        ## The covariate files were written by the catalog arm's Module 01 run,
        ## so they carry ITS prefix, not the group's. locus_io.R needs this
        ## explicitly: its fallback ternary would look for TOPMed_LIBD.AA for an
        ## all_individuals.AA cell, and that file does not exist here.
        covar_prefix     = arm$covar_prefix,
        n_donors_pooled  = nrow(pooled_donors),
        n_donors_cell    = nrow(cell_donors),
        ## Also written under the Module 01 field name. Module 02's Stage 00
        ## reads n_donors and donor_checksum straight out of its upstream
        ## manifest, and an 01b run has to be a drop-in substitute there.
        n_donors         = nrow(cell_donors),
        donor_checksum   = donor_checksum(cell_donors$FID),
        cis_window_bp    = load_config("thresholds")$cis$window_bp,
        genotype_pcs     = paste(
            load_config("covariates")$estimation_cells$genotype_pcs,
            collapse = ",")
    ))

## Same layout as a Module 01 run, because 02 and 03 read that layout.
vmr_dir <- file.path(run$dir, "vmr")
dir.create(vmr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run$dir, "covs"), recursive = TRUE, showWarnings = FALSE)

write_atomic(cell_donors, file.path(vmr_dir, "donors_plink.txt"),
             col.names = FALSE)

## The locus set is copied verbatim, never recomputed. Both cells on a catalog
## get byte-identical copies; stage 04 checksums them against the source.
##
## vmr.bed sizes the extraction array; vmr_catalog.tsv is what Module 02 reads
## for vmr_id and vmr_set_id. Copying both is what makes an 01b run a drop-in
## substitute for a Module 01 run at the 02/03 call site.
for (f in c("vmr.bed", "vmr_catalog.tsv")) {
    src <- file.path(vmr_run_dir, "vmr", f)
    if (!file.exists(src)) stop("Module 01 run lacks vmr/", f, ": ", vmr_run_dir)
    file.copy(src, file.path(vmr_dir, f), overwrite = FALSE)
    if (!file.exists(file.path(vmr_dir, f))) {
        stop("Failed to copy ", f, " from ", vmr_run_dir)
    }
}

message("[cell] ", cell, "/", region, ": ", nrow(cell_donors), " of ",
        nrow(pooled_donors), " pooled donors")

## The launcher reads this line.
cat(run$run_id, "\n", sep = "")
