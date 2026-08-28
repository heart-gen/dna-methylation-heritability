#### 01_vmr_catalog / 04_turnover: old-vs-new comparison and QC tables ####
##
## AGENTS.md 7.1 requires an old-versus-new turnover table, an array-coverage
## comparison, and a technical QC / exclusion table.
##
## AGENTS.md 8: the legacy catalog is a comparison baseline ONLY. It is invalid
## for scientific use (V1: caudate chr1 had 153/153 donor rows misaligned), so
## nothing here licenses reusing a legacy number. The point of the comparison is
## to quantify how much moved, and to catch a change so large it indicates a new
## bug rather than the known repair.
##
## Usage:
##   Rscript 04_turnover.R --cohort AA --region caudate --run-id ID

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
})

opts <- parse_v2_args(require = c("cohort", "region", "run_id"))
cohort <- opts$cohort; region <- opts$region

th <- load_config("thresholds")
module_root <- file.path(V2_ROOT, "01_vmr_catalog")
run_dir <- file.path(module_root, "_m", "runs", opts$run_id)
vmr_dir <- file.path(run_dir, "vmr")
qc_dir  <- file.path(run_dir, "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

new_vmr <- fread(file.path(vmr_dir, "vmr_catalog.tsv"))
vmr_set_id <- new_vmr$vmr_set_id[1]

## The legacy catalog for this arm x region, if it is still on disk.
legacy_path <- if (cohort == "AA") {
    file.path(V2_ROOT, "vmr-analysis", region, "_m", "vmr.bed")
} else {
    file.path(V2_ROOT, "vmr-analysis", "all_individuals", region, "_m", "vmr.bed")
}

to_gr <- function(dt) GRanges(dt$chr, IRanges(dt$start, dt$end))

## ------------------------------------------------------------ turnover table

if (file.exists(legacy_path)) {
    legacy <- fread(legacy_path, header = FALSE, colClasses = list(character = 1),
                    col.names = c("chr", "start", "end"))
    ## The legacy catalogs write bare chromosome names ("1", "X"); v2 writes
    ## "chr1". Normalize before comparing, or every region looks novel.
    legacy[, chr := paste0("chr", sub("^chr", "", chr))]

    ## Legacy catalogs include unmasked sex chromosomes (V4). Compare on the
    ## primary autosomes only, or the turnover number is dominated by a policy
    ## change rather than by the V1 repair.
    legacy_auto <- legacy[chr %in% paste0("chr", chrom_order(include_sex = FALSE))]
    n_legacy_sex <- nrow(legacy) - nrow(legacy_auto)

    ## A smoke run covers only some chromosomes. Restrict the legacy side to the
    ## chromosomes v2 actually produced, so the comparison is like-for-like
    ## rather than reporting every unbuilt chromosome as "lost".
    v2_chroms <- unique(new_vmr$chr)
    legacy_auto <- legacy_auto[chr %in% v2_chroms]
    if (length(v2_chroms) < length(chrom_order(include_sex = FALSE))) {
        message("[turnover] v2 catalog covers ", length(v2_chroms),
                " chromosome(s); comparing on those only.")
    }

    g_new <- to_gr(new_vmr); g_old <- to_gr(legacy_auto)
    ov <- findOverlaps(g_new, g_old)

    inter <- sum(width(GenomicRanges::intersect(g_new, g_old, ignore.strand = TRUE)))
    union_w <- sum(width(GenomicRanges::union(g_new, g_old, ignore.strand = TRUE)))

    turnover <- data.table(
        cohort = cohort, region = region, vmr_set_id = vmr_set_id,
        chromosomes_compared = paste(sort(unique(new_vmr$chr)), collapse = ","),
        n_v2 = nrow(new_vmr),
        n_legacy_autosomal = nrow(legacy_auto),
        n_legacy_sex_chrom_dropped = n_legacy_sex,
        n_v2_overlapping_legacy = length(unique(queryHits(ov))),
        n_v2_novel = nrow(new_vmr) - length(unique(queryHits(ov))),
        n_legacy_lost = nrow(legacy_auto) - length(unique(subjectHits(ov))),
        jaccard_bp = if (union_w > 0) inter / union_w else NA_real_,
        frac_v2_novel = 1 - length(unique(queryHits(ov))) / nrow(new_vmr),
        legacy_path = sub(paste0("^", V2_ROOT, "/"), "", legacy_path),
        legacy_status = "legacy_invalid_for_prediction_accuracy")

    print(turnover)

    gate <- th$gates$max_vmr_turnover
    if (!is.null(gate) && !is.na(gate)) {
        if (turnover$frac_v2_novel > gate) {
            stop("VMR turnover ", round(turnover$frac_v2_novel, 3),
                 " exceeds the locked threshold ", gate,
                 ".\n  AGENTS.md 14: stop and request PI direction.")
        }
    } else {
        message("[gate] max_vmr_turnover is not locked; reporting turnover ",
                "without gating (AGENTS.md 12).")
    }
} else {
    turnover <- data.table(
        cohort = cohort, region = region, vmr_set_id = vmr_set_id,
        n_v2 = nrow(new_vmr), legacy_path = legacy_path,
        legacy_status = "legacy_catalog_not_on_disk")
    message("[turnover] no legacy catalog at ", legacy_path,
            "; recording v2 counts only")
}
write_atomic(turnover, file.path(qc_dir, "vmr_turnover.tsv"))

## ------------------------------------------------------ array-coverage table
##
## One of the manuscript's stated contributions is WGBS coverage of methylation
## outside array-accessible CpGs (AGENTS.md 2.2). This quantifies it.
##
## The array universes are reference annotations describing a platform we did
## not run -- coordinates only, in the same category as repeat-masker or the
## ENCODE blacklist. No array intensity, sample, or methylation value enters
## this computation, and this stage runs after VMR calling and writes only to
## qc/, so no VMR is defined or filtered using array positions.
##
## Two platforms are reported. 450K is the required main-figure comparator;
## EPIC is the stricter supplemental one (~866k vs ~486k probes), so the
## off-array fraction against EPIC must be the smaller of the two.
##
## Schema of both universes (built by inputs/supportfiles/_h/01_build_array_universe.R):
##   probe_id  chrom  pos_1based

membership <- fread(file.path(vmr_dir, "cpg_vmr_membership.tsv"))

support_dir <- file.path(V2_ROOT, "inputs", "supportfiles", "_m")
ARRAY_UNIVERSES <- list(
    list(platform = "450K", file = "450k_universe_hg38.tsv.gz", required = TRUE),
    list(platform = "EPIC", file = "epic_universe_hg38.tsv.gz", required = FALSE)
)

g_cpg <- GRanges(membership$chr, IRanges(membership$cpg_pos, width = 1))
g_vmr <- to_gr(new_vmr)

coverage <- rbindlist(lapply(ARRAY_UNIVERSES, function(spec) {
    path <- file.path(support_dir, spec$file)

    ## A silent skip here is what let this panel go missing across six accepted
    ## runs: the old code wrote a placeholder note and exited 0. The required
    ## comparator now fails the stage instead.
    if (!file.exists(path)) {
        if (spec$required) {
            stop("Required array universe not found: ", path,
                 "\nAGENTS.md 2.2 and 11 (Figure 1) need the off-array ",
                 "coverage comparison. Build it with:\n",
                 "  Rscript inputs/supportfiles/_h/01_build_array_universe.R ",
                 "--platform ", spec$platform)
        }
        message("[coverage] optional universe absent, skipping: ", spec$platform)
        return(NULL)
    }

    arr <- fread(path)
    stopifnot(all(c("chrom", "pos_1based") %in% names(arr)))
    g_arr <- GRanges(arr$chrom, IRanges(arr$pos_1based, width = 1))

    cpg_on_array <- countOverlaps(g_cpg, g_arr) > 0
    vmr_any_array <- countOverlaps(g_vmr, g_arr) > 0

    data.table(
        cohort = cohort, region = region, vmr_set_id = vmr_set_id,
        array_platform = spec$platform,
        n_array_probes = nrow(arr),
        n_vmr_cpgs = nrow(membership),
        n_vmr_cpgs_on_array = sum(cpg_on_array),
        frac_vmr_cpgs_off_array = mean(!cpg_on_array),
        n_vmrs = nrow(new_vmr),
        n_vmrs_with_no_array_cpg = sum(!vmr_any_array),
        frac_vmrs_invisible_to_array = mean(!vmr_any_array),
        array_manifest = basename(path))
}), fill = TRUE)

## Sanity: more probes must not yield more off-array CpGs.
if (all(c("450K", "EPIC") %in% coverage$array_platform)) {
    f450 <- coverage[array_platform == "450K", frac_vmr_cpgs_off_array]
    fepic <- coverage[array_platform == "EPIC", frac_vmr_cpgs_off_array]
    if (fepic > f450) {
        stop("EPIC off-array fraction (", signif(fepic, 4), ") exceeds 450K (",
             signif(f450, 4), "); EPIC is a superset platform, so this ",
             "indicates a coordinate or build mismatch.")
    }
}

write_atomic(coverage, file.path(qc_dir, "array_coverage.tsv"))

## ------------------------------------------- technical QC and exclusion table

prep_files <- list.files(run_dir, pattern = "^prepare_summary\\.tsv$",
                         recursive = TRUE, full.names = TRUE)
qc_rows <- rbindlist(lapply(prep_files, function(f) {
    d <- fread(f)
    out <- as.list(d$value); names(out) <- d$field
    as.data.table(out)
}), fill = TRUE)

## Sex chromosomes prepared but held out of the primary catalog (V4).
excluded_dir <- file.path(run_dir, "excluded")
excl <- if (dir.exists(excluded_dir)) {
    sex_prep <- list.files(excluded_dir, pattern = "^prepare_summary\\.tsv$",
                           recursive = TRUE, full.names = TRUE)
    rbindlist(lapply(sex_prep, function(f) {
        d <- fread(f); out <- as.list(d$value); names(out) <- d$field
        cbind(as.data.table(out),
              exclusion_reason = "sex_chromosome_no_ct_mask_not_in_primary_catalog")
    }), fill = TRUE)
} else data.table()

write_atomic(qc_rows, file.path(qc_dir, "technical_qc.tsv"))
if (nrow(excl) > 0) write_atomic(excl, file.path(qc_dir, "exclusions.tsv"))

message("[done] QC tables written to ", qc_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
