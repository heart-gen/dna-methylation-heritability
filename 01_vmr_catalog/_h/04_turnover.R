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
source(file.path(V2_ROOT, "01_vmr_catalog", "_h", "exclusion_ledger.R"))

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
##
## This glob deliberately reaches into excluded/ as well: `is_primary_chrom`
## distinguishes those rows and Figure 1 filters on it, so a sex chromosome
## cannot inflate the assayed-CpG count while still being visible here.
prep_files <- list.files(run_dir, pattern = "^prepare_summary\\.tsv$",
                         recursive = TRUE, full.names = TRUE)
read_field_table <- function(f) {
    d <- fread(f)
    out <- as.list(d$value); names(out) <- d$field
    as.data.table(out)
}
qc_rows <- rbindlist(lapply(prep_files, read_field_table), fill = TRUE)
write_atomic(qc_rows, file.path(qc_dir, "technical_qc.tsv"))

## ============================================================================
## F14: the exclusion table, the chromosome-policy manifest, and the accounting
## ============================================================================
##
## AGENTS.md 7.1 requires a "technical QC and exclusion table" and, for the
## chromosome policy, that sex chromosomes are "reported separately or excluded
## with an explicit manifest". Before this block, `qc/exclusions.tsv` was written
## only when a sex-chromosome preparation happened to be present -- and
## step_1x.sh is opt-in (`WITH_SEX=1`), so all six accepted runs have an empty
## excluded/ and no exclusions.tsv at all. Worse, even when present it would have
## held only sex chromosomes: no donor, CpG or candidate-VMR exclusion was
## recorded anywhere but a SLURM log.
##
## Nothing here decides an exclusion. The upstream stages now write their own
## decisions as ledger parts; this assembles them, and balances them against the
## surviving counts so a denominator cannot be quietly wrong (AGENTS.md 11, 14).

## ------------------------------------------ chromosome policy manifest
##
## Written unconditionally, and derived from config rather than from what happens
## to be on disk. That is the difference that matters: the manifest states that
## X and Y are held out of the primary catalog and why, whether or not anyone
## chose to run step_1x.sh for this run.
primary_chroms <- as.character(chrom_order(include_sex = FALSE))
sex_chrom_ids <- sex_chroms(th)

status_file <- file.path(vmr_dir, "chromosome_status.tsv")
chrom_status <- if (file.exists(status_file)) {
    fread(status_file, colClasses = list(character = "chrom"))
} else {
    data.table(chrom = character(), status = character())
}

cutoffs_file <- file.path(vmr_dir, "sd_cutoffs.tsv")
cutoffs <- if (file.exists(cutoffs_file)) fread(cutoffs_file) else data.table()

chrom_prepared <- function(cc) {
    any(dir.exists(c(
        file.path(run_dir, "cpg", paste0("chr_", cc)),
        file.path(run_dir, "excluded", "cpg", paste0("chr_", cc)))))
}

chrom_manifest <- rbindlist(lapply(c(primary_chroms, sex_chrom_ids), function(cc) {
    in_policy <- cc %in% primary_chroms
    st <- chrom_status$status[match(cc, chrom_status$chrom)]
    prepared <- chrom_prepared(cc)
    disposition <- if (!in_policy) {
        "excluded_from_primary_catalog_by_chromosome_policy"
    } else if (!is.na(st) && st == "completed") {
        "in_primary_catalog"
    } else if (!is.na(st)) {
        st
    } else if (!prepared) {
        "not_prepared_in_this_run"
    } else {
        "prepared_but_catalog_status_unrecorded"
    }
    reason <- if (!in_policy) {
        ## The substantive reason, not the policy restated: no C->T SNP mask
        ## exists for X/Y, so those CpGs are unmasked and the legacy caudate
        ## catalog carried a 3x excess of sex-chromosome VMRs (V4).
        "no_ct_snp_mask_available_unmasked_cpgs_v4"
    } else if (disposition == "in_primary_catalog") {
        NA_character_
    } else {
        disposition
    }
    data.table(
        cohort = cohort, region = region, vmr_set_id = vmr_set_id,
        chrom = cc, in_primary_policy = in_policy,
        ct_mask_available = has_ct_mask(cc),
        prepared_in_this_run = prepared,
        disposition = disposition, reason = reason,
        n_vmrs_in_primary_catalog = if (nrow(cutoffs) > 0) {
            as.integer(cutoffs$n_vmrs[match(paste0("chr", cc), cutoffs$chr)])
        } else NA_integer_,
        reported_separately_under = if (in_policy) NA_character_ else {
            if (prepared) "_m/runs/{RUN_ID}/excluded/" else "not_run"
        })
}), use.names = TRUE)
write_atomic(chrom_manifest, file.path(qc_dir, "chromosome_policy_manifest.tsv"))

## ------------------------------------------------- assemble qc/exclusions.tsv
ledger_files <- list.files(run_dir, pattern = "^exclusions_[a-z_]+\\.tsv$",
                           recursive = TRUE, full.names = TRUE)
ledger_parts <- lapply(ledger_files, function(f) {
    d <- fread(f, colClasses = list(character = "chrom"))
    if (nrow(d) == 0) return(NULL)
    rel <- sub(paste0("^", run_dir, "/"), "", f)
    d[, catalog_tree := if (startsWith(rel, "excluded/")) "excluded" else "primary"]
    d[, source_file := rel][]
})
ledger <- rbindlist(c(list(excl_empty()[, `:=`(catalog_tree = character(),
                                               source_file = character())]),
                      ledger_parts), use.names = TRUE, fill = TRUE)

## Donor and CpG parts are written per chromosome. Donor eligibility does not
## vary by chromosome, so 22 identical rows per donor are collapsed to one --
## and if they ever are not identical, the rows stay separate and say so.
donor_rows <- ledger[unit_type %in% c("donor", "phenotype_row") &
                     catalog_tree == "primary"]
other_rows <- ledger[!(unit_type %in% c("donor", "phenotype_row") &
                       catalog_tree == "primary")]
if (nrow(donor_rows) > 0) {
    ## source_file names the chromosome directory, which would defeat the
    ## collapse. Glob it, so the provenance survives as "which stage part" while
    ## staying chromosome-independent -- the same property the rows themselves
    ## have.
    donor_rows[, source_file := sub("chr_[0-9XY]+/", "chr_*/", source_file)]
    donor_rows <- collapse_uniform_chrom(
        donor_rows, primary_chroms, by_extra = c("catalog_tree", "source_file"))
}

## The chromosome policy joins the ledger as chromosome-level rows, so a single
## table answers "what did not make it into this catalog".
chrom_rows <- chrom_manifest[!is.na(reason), excl_rows(
    stage = "chromosome_policy", unit_type = "chromosome",
    exclusion_reason = reason, unit_id = paste0("chr", chrom), chrom = chrom)]
if (nrow(chrom_rows) > 0) {
    chrom_rows[, `:=`(catalog_tree = "primary",
                      source_file = "config/thresholds.yml:chromosomes")]
}

exclusions <- rbindlist(list(donor_rows, other_rows, chrom_rows),
                        use.names = TRUE, fill = TRUE)
exclusions[, `:=`(cohort = cohort, region = region, vmr_set_id = vmr_set_id,
                  run_id = opts$run_id)]
setcolorder(exclusions, c("cohort", "region", "run_id", "vmr_set_id",
                          EXCLUSION_COLS, "catalog_tree", "source_file"))
write_atomic(exclusions, file.path(qc_dir, "exclusions.tsv"))

## -------------------------------------------------- entered = survived + excluded
##
## The accounting table is the claim; the ledger is the evidence. A row whose
## inputs are present and whose arithmetic does not close stops the stage, since
## AGENTS.md 14 makes a degenerate or undocumented denominator a stop condition.
## A row whose inputs are absent is recorded as unbalanced with a note, never as
## balanced -- that is the case for a QC-refresh run whose source predates these
## ledger parts.
num <- function(dt, field_name) {
    if (is.null(dt) || nrow(dt) == 0 || !field_name %in% names(dt)) {
        return(NA_integer_)
    }
    v <- suppressWarnings(as.integer(dt[[field_name]]))
    if (length(v) == 0) NA_integer_ else v
}
one <- function(x, what) {
    u <- unique(x[!is.na(x)])
    if (length(u) == 0) return(NA_integer_)
    if (length(u) > 1) {
        stop(what, " differs across autosomes (", paste(u, collapse = ", "),
             "). Donor eligibility must not depend on the chromosome.")
    }
    u
}

prim <- if ("is_primary_chrom" %in% names(qc_rows)) {
    qc_rows[as.logical(is_primary_chrom) %in% TRUE]
} else qc_rows

summ_file <- file.path(vmr_dir, "summarize_summary.tsv")
summ <- if (file.exists(summ_file)) {
    d <- fread(summ_file); out <- as.list(d$value); names(out) <- d$field
    as.data.table(out)
} else data.table()

analyze_files <- list.files(file.path(run_dir, "pca"),
                            pattern = "^analyze_summary\\.tsv$",
                            recursive = TRUE, full.names = TRUE)
analyze <- if (length(analyze_files) > 0) {
    rbindlist(lapply(analyze_files, read_field_table), fill = TRUE)
} else data.table()

no_ledger <- paste("upstream ledger fields absent; rerun the module to populate",
                   "them (a QC-refresh run copies a pre-F14 source)")

acc <- list()

## Phenotype rows: the whole table, split into this region's candidates and the
## rows that belong to another region's run.
n_cand <- one(num(prim, "n_donor_candidate_rows"), "n_donor_candidate_rows")
n_other <- one(num(prim, "n_phenotype_rows_other_region"),
               "n_phenotype_rows_other_region")
acc[["pheno"]] <- accounting_row(
    "phenotype_row", "phenotype_table_all_regions",
    if (is.na(n_cand) || is.na(n_other)) NA_integer_ else n_cand + n_other,
    n_cand, n_other,
    note = if (is.na(n_cand)) no_ledger else
        "survivors are this region's candidate donors, not catalog donors")

## Donors: candidates for this region against the catalog's donor set.
n_donors_catalog <- num(summ, "n_donors")
acc[["donor"]] <- accounting_row(
    "donor", "autosomes_all", n_cand, n_donors_catalog,
    one(num(prim, "n_donors_excluded"), "n_donors_excluded"),
    note = if (is.na(n_cand)) no_ledger else
        if (is.na(n_donors_catalog)) "vmr/summarize_summary.tsv absent" else
            NA_character_)

## Donors again, at the residualization step: prepared against residualized.
if (nrow(analyze) > 0 && "n_donors_prepared" %in% names(analyze)) {
    acc[["donor_resid"]] <- accounting_row(
        "donor", "residualization_autosomes_total",
        sum(num(analyze, "n_donors_prepared")),
        sum(num(analyze, "n_donors")),
        sum(num(analyze, "n_donors_dropped_missing_snp_pcs")),
        note = "summed over chromosomes; a nonzero excluded count means the VMR cutoff and the VMR phenotypes used different donor sets")
}

## CpGs, per chromosome and in total.
if ("n_cpgs_input" %in% names(prim)) {
    cpg_in <- num(prim, "n_cpgs_input"); cpg_out <- num(prim, "n_cpgs")
    cpg_ct <- num(prim, "n_cpgs_excluded_ct_snp")
    cpg_lc <- num(prim, "n_cpgs_excluded_low_coverage")
    acc[["cpg_chrom"]] <- rbindlist(lapply(seq_len(nrow(prim)), function(i) {
        accounting_row("cpg", paste0("chr", prim$chrom[i]), cpg_in[i],
                       cpg_out[i], cpg_ct[i] + cpg_lc[i])
    }))
    acc[["cpg_total"]] <- accounting_row(
        "cpg", "autosomes_total", sum(cpg_in), sum(cpg_out),
        sum(cpg_ct) + sum(cpg_lc))
} else {
    acc[["cpg_total"]] <- accounting_row("cpg", "autosomes_total", NA_integer_,
                                         sum(num(prim, "n_cpgs")), NA_integer_,
                                         note = no_ledger)
}

## Candidate VMRs: what regionFinder3() proposed against what min_cpgs kept.
acc[["vmr"]] <- accounting_row(
    "vmr_candidate", "autosomes_all", num(summ, "n_vmr_candidates"),
    nrow(new_vmr), num(summ, "n_vmr_candidates_excluded"),
    note = if (is.na(num(summ, "n_vmr_candidates"))) no_ledger else NA_character_)

## Chromosomes: every chromosome in the genome policy is in exactly one bucket.
acc[["chrom"]] <- accounting_row(
    "chromosome", "genome_policy", nrow(chrom_manifest),
    sum(chrom_manifest$disposition == "in_primary_catalog"),
    sum(chrom_manifest$disposition != "in_primary_catalog"))

accounting <- rbindlist(acc, use.names = TRUE, fill = TRUE)
accounting[, `:=`(cohort = cohort, region = region, run_id = opts$run_id,
                  vmr_set_id = vmr_set_id)]
setcolorder(accounting, c("cohort", "region", "run_id", "vmr_set_id",
                          ACCOUNTING_COLS))
write_atomic(accounting, file.path(qc_dir, "exclusion_accounting.tsv"))
print(accounting[, .(unit_type, scope, n_entered, n_survived, n_excluded,
                     balanced)])

broken <- accounting[balanced == FALSE & !is.na(n_entered) &
                     !is.na(n_survived) & !is.na(n_excluded)]
if (nrow(broken) > 0) {
    stop("Exclusion accounting does not balance for ", nrow(broken), " row(s): ",
         paste(broken$unit_type, broken$scope, sep = "/", collapse = ", "),
         "\n  entered must equal survived + excluded. AGENTS.md 14 lists an ",
         "undocumented denominator as a stop condition.")
}
unrecorded <- accounting[balanced == FALSE]
if (nrow(unrecorded) > 0) {
    message("[exclusions] ", nrow(unrecorded), " accounting row(s) could not be ",
            "balanced because their inputs are absent: ",
            paste(unrecorded$unit_type, unrecorded$scope, sep = "/",
                  collapse = ", "),
            "\n  This is expected for a QC refresh of a pre-F14 source run and ",
            "is recorded in the table, not swallowed.")
}

message("[done] QC tables written to ", qc_dir, ": technical_qc.tsv, ",
        "exclusions.tsv (", nrow(exclusions), " rows), ",
        "exclusion_accounting.tsv, chromosome_policy_manifest.tsv")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
