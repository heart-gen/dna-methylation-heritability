#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- decide which CpGs are tested, and record why not ####
##
## Usage:
##   Rscript _h/01_prepare_cpg_set.R --run-id cmb-AA-caudate-20260823
##
## This stage owns the DENOMINATOR. Every burden fraction downstream divides by
## `n_tested_cpgs` from the table this script writes, never by the number of
## CpGs a VMR contains. AGENTS.md 7.5 requires tested CpGs to be reported apart
## from prepared-but-untested ones and every concordance denominator audited,
## and REVISION_GUIDE R8 records degenerate denominators as a real defect of the
## legacy module.
##
## Module 01's layout is per chromosome: `cpg/chr_{N}/cpg_meth.phen` is a
## donors x CpGs matrix whose column names are bare positions, and
## `vmr/cpg_vmr_membership.tsv` keys CpGs by (chr, cpg_pos) with no ID column.
## The stable identifier is therefore constructed here, once, as `chr:pos`, and
## everything downstream uses it.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "05_cpg_meqtl_burden"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mval("cohort"); region <- mval("region")
smoke <- identical(mval("smoke_run"), "TRUE")
meqtl <- load_config("meqtl_parameters")
thresholds <- load_config("thresholds")

cat_dir <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs",
                     mval("upstream_vmr_catalog_run_id"))
if (!dir.exists(cat_dir)) stop("Upstream 01 run not found: ", cat_dir)

## ------------------------------------------------------- corrected membership
membership <- fread(file.path(cat_dir, "vmr", "cpg_vmr_membership.tsv"))
required <- c("chr", "cpg_pos", "vmr_id")
if (!all(required %in% names(membership))) {
    stop("Membership table is missing required columns; expected ",
         paste(required, collapse = ", "), "; found ",
         paste(names(membership), collapse = ", "))
}
membership[, cpg_id := paste(chr, cpg_pos, sep = ":")]

## The catalog table carries the vmr_set_id; the membership table does not.
catalog <- fread(file.path(cat_dir, "vmr", "vmr_catalog.tsv"))
if (!identical(as.character(catalog$vmr_set_id[1]), mval("vmr_set_id"))) {
    stop("Module 01 catalog cites vmr_set_id ", catalog$vmr_set_id[1],
         " but this run was opened against ", mval("vmr_set_id"))
}

## Restrict to the autosomes the chromosome policy declares. A CpG on a
## chromosome the policy excludes is not a QC failure -- it was never in scope.
primary <- paste0("chr", config_get(thresholds, "chromosomes.primary"))
off_policy <- membership[!chr %in% primary]
membership <- membership[chr %in% primary]
## A smoke run may cover a subset of the autosomes. Record that subset in the
## manifest: stage 02b reconciles the per-chromosome array against the config's
## full autosome list, and without a recorded restriction it would -- correctly
## -- refuse the run for 21 unaccounted chromosomes. The restriction is only
## ever honoured for `smoke_run = TRUE`, so a production run cannot narrow its
## own denominator this way.
smoke_chroms <- NULL
if (smoke && !is.null(opts$smoke_chroms)) {
    keep_chr <- paste0("chr", trimws(strsplit(opts$smoke_chroms, ",")[[1]]))
    membership <- membership[chr %in% keep_chr]
    smoke_chroms <- paste(keep_chr, collapse = ",")
    message("[05] SMOKE: restricted to ", paste(keep_chr, collapse = ","))
}
if (nrow(membership) == 0) stop("No member CpGs on the in-policy chromosomes")

## --------------------------------------------------------------- exclusions
## Every exclusion is recorded with a reason, and the reasons are mutually
## exclusive and ordered, so a CpG appears exactly once in the audit.
dir.create(file.path(run_dir, "results", "tested_meth"),
           recursive = TRUE, showWarnings = FALSE)

excl <- data.table(cpg_id = character(), reason = character())
add_excl <- function(ids, reason) {
    ids <- setdiff(unique(ids), excl$cpg_id)
    if (length(ids) > 0) {
        excl <<- rbind(excl, data.table(cpg_id = ids, reason = reason))
    }
}

qc <- meqtl$cpg_qc
donors_seen <- NULL
present <- character()

## Non-residualized CpG methylation is the primary source (meqtl_parameters.yml).
## Residualized values are permitted only as a negative-control sensitivity:
## residualizing on PCs derived from the same CpGs removes part of the genetic
## signal being tested.
for (cc in sort(unique(membership$chr))) {
    chrom_n <- sub("^chr", "", cc)
    meth_f <- file.path(cat_dir, "cpg", paste0("chr_", chrom_n), "cpg_meth.phen")
    if (!file.exists(meth_f)) {
        add_excl(membership[chr == cc, cpg_id], "chromosome_matrix_missing")
        next
    }
    ## header = TRUE is not optional. Every CpG column name is a bare position,
    ## so fread's autodetection sees an all-numeric first row and decides the
    ## file has no header -- which silently turns the literal string "FID" into
    ## a 154th donor and shifts every column name by one.
    ## Read the HEADER first and select only the member columns. These matrices
    ## are donors x every CpG on the chromosome -- chr1 is 5.2 GB and roughly
    ## two million columns -- and reading one whole is not merely slow: data
    ## .table walks its columns with vapply during the read, so a two-million
    ## column data.table segfaults the R session ("recursive gc invocation").
    ## Only a few thousand columns per chromosome are ever VMR members, so the
    ## narrowed read is both the fast path and the only one that survives.
    hdr <- names(fread(meth_f, header = TRUE, nrows = 0L))
    col_ids_all <- paste(cc, hdr[-(1:2)], sep = ":")
    mem_ids <- membership[chr == cc, cpg_id]
    add_excl(setdiff(mem_ids, col_ids_all), "not_in_methylation_matrix")

    shared <- intersect(mem_ids, col_ids_all)
    if (length(shared) == 0) next
    sel <- match(shared, col_ids_all) + 2L

    ## fread cannot open this file at all -- not even with `select`. At roughly
    ## two million columns it segfaults ("memory not mapped") while setting up
    ## per-column state, before any row is parsed, at a peak RSS of ~6 GB. The
    ## header-only read above survives because it stops at one line.
    ##
    ## So the column subset is taken by `cut`, which streams and does not care
    ## how wide the line is, and only the resulting narrow file is handed to
    ## fread. The field list is a few thousand indices, far short of ARG_MAX.
    narrow_f <- tempfile(pattern = paste0("cpgsel_", cc, "_"), fileext = ".tsv")
    on.exit(unlink(narrow_f), add = TRUE)
    fields <- paste(c(1L, 2L, sort(sel)), collapse = ",")
    cut_status <- system2("cut", c("-f", shQuote(fields), shQuote(meth_f)),
                          stdout = narrow_f)
    if (!identical(cut_status, 0L)) {
        stop("cut failed (status ", cut_status, ") extracting ", length(sel),
             " columns from ", meth_f)
    }
    meth <- fread(narrow_f, header = TRUE, colClasses = list(character = 1:2))
    unlink(narrow_f)
    donors <- as.character(meth[[1]])          # FID is the donor ID
    assert_no_dups(donors, paste0("donors in ", meth_f))
    if (is.null(donors_seen)) {
        donors_seen <- donors
        assert_expected_n(length(donors), cohort, region)
    } else if (!identical(donors, donors_seen)) {
        stop("Donor set differs between chromosomes: ", cc,
             " does not match the first chromosome read.")
    }

    ## `select` returns columns in file order, so re-key from the names fread
    ## actually gave back rather than assuming the requested order was kept.
    got_ids <- paste(cc, names(meth)[-(1:2)], sep = ":")
    mat <- as.matrix(meth[, match(shared, got_ids) + 2L, with = FALSE])

    frac_nonmissing <- colMeans(!is.na(mat))
    add_excl(shared[frac_nonmissing < qc$min_fraction_samples_passing_coverage],
             "below_min_fraction_samples_passing_coverage")

    if (isTRUE(qc$require_nonzero_variance)) {
        v <- matrixStats::colVars(mat, na.rm = TRUE)
        add_excl(shared[!is.finite(v) | v == 0], "zero_variance")
    }
    present <- c(present, shared)

    ## Emit the tested methylation for this chromosome, donors x CpGs, keyed by
    ## the chr:pos identifier. 01b then builds the tensorqtl BED from a few
    ## thousand columns instead of re-reading a ~2M-column matrix in Python --
    ## which gets OOM-killed. data.table has already done the wide read here;
    ## doing it twice, in two languages, is both slow and a place for the two
    ## CpG-identifier conventions to drift apart.
    keep_now <- setdiff(shared, excl$cpg_id)
    if (length(keep_now)) {
        tested_mat <- cbind(
            data.table(FID = donors),
            as.data.table(mat[, match(keep_now, shared), drop = FALSE]))
        setnames(tested_mat, c("FID", keep_now))
        ## Plain .tsv, not .tsv.gz: write_atomic() writes to a "<path>.tmp.<pid>"
        ## file and renames, so fwrite never sees a .gz extension and never
        ## compresses. A .gz name on uncompressed bytes is worse than no
        ## compression -- the Python reader fails on the magic number.
        write_atomic(tested_mat,
                     file.path(run_dir, "results", "tested_meth",
                               paste0(cc, ".tsv")))
    }
    rm(meth, mat); invisible(gc(FALSE))
}
if (is.null(donors_seen)) stop("No chromosome methylation matrix could be read")

## ENCODE blacklist. The config has declared this since the module was written
## and the previous draft silently ignored it -- a declared QC rule that never
## ran is worse than one that was never declared.
if (isTRUE(qc$exclude_blacklist)) {
    bl_f <- file.path(repo_root(), "inputs", "supportfiles", "_m",
                      "hg38-blacklist.v2.bed.gz")
    if (!file.exists(bl_f)) stop("exclude_blacklist is TRUE but no blacklist at ", bl_f)
    bl <- fread(cmd = paste("zcat", shQuote(bl_f)), header = FALSE, select = 1:3,
                col.names = c("chr", "start", "end"))
    setkey(bl, chr, start, end)
    hit <- foverlaps(membership[, .(chr, start = cpg_pos, end = cpg_pos, cpg_id)],
                     bl, type = "within", nomatch = NULL)
    add_excl(hit$cpg_id, "encode_blacklist")
}

## C>T SNPs at the CpG: a common C>T destroys the CpG and produces an apparent
## methylation difference that is really a genotype readout, which would
## manufacture meQTLs.
if (isTRUE(qc$exclude_common_ct_snp_at_cpg)) {
    ## One file per chromosome, each a bare list of positions. Address them
    ## through the config template rather than by listing the directory: that
    ## directory also holds the job scripts and log directories that produced
    ## the lists, and globbing it reads those as data.
    ##
    ## Masks exist for autosomes only (00_shared/chrom.R::has_ct_mask), so a
    ## chromosome without one is recorded as unmasked rather than treated as
    ## having no C>T SNPs.
    unmasked <- character()
    for (cc in sort(unique(membership$chr))) {
        chrom_n <- sub("^chr", "", cc)
        if (!has_ct_mask(chrom_n)) { unmasked <- c(unmasked, cc); next }
        ct_f <- resolve_path("ct_snp_template", chrom = chrom_n)
        if (!file.exists(ct_f)) { unmasked <- c(unmasked, cc); next }
        ct_pos <- fread(ct_f, header = FALSE)[[1]]
        ## Match on position, not on ID: the C>T lists and the BSobj use
        ## different identifier conventions and an ID join matches nothing.
        add_excl(intersect(membership[chr == cc, cpg_id],
                           paste(cc, ct_pos, sep = ":")),
                 "common_ct_snp_at_cpg")
    }
    if (length(unmasked)) {
        write_atomic(data.table(chrom = unmasked, reason = "no_ct_mask_available"),
                     file.path(run_dir, "results", "unmasked-chromosomes.tsv"))
        warning("No C>T mask for: ", paste(unmasked, collapse = ", "),
                ". Recorded in results/unmasked-chromosomes.tsv.", call. = FALSE)
    }
}

tested <- setdiff(present, excl$cpg_id)
if (length(tested) == 0) stop("No CpGs survived QC; check the exclusion reasons")

## ------------------------------------------------------------------ outputs
out <- file.path(run_dir, "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

write_atomic(membership[cpg_id %in% tested], file.path(out, "tested-cpg-membership.tsv"))
write_atomic(excl, file.path(out, "excluded-cpgs.tsv"))
if (nrow(off_policy)) {
    write_atomic(off_policy[, .N, by = chr],
                 file.path(out, "off-policy-chromosome-cpgs.tsv"))
}

## The denominator table. Every downstream burden fraction must use
## n_tested_cpgs from here, never n_member_cpgs.
denom <- membership[, .(n_member_cpgs = .N), by = vmr_id]
denom <- merge(denom,
               membership[cpg_id %in% tested, .(n_tested_cpgs = .N), by = vmr_id],
               by = "vmr_id", all.x = TRUE)
denom[is.na(n_tested_cpgs), n_tested_cpgs := 0L]
write_atomic(denom, file.path(out, "vmr-cpg-denominators.tsv"))

summary_dt <- data.table(
    region = region, population = cohort, vmr_set_id = mval("vmr_set_id"),
    n_member_cpgs = nrow(membership),
    n_tested_cpgs = length(tested),
    n_excluded_cpgs = nrow(excl),
    n_off_policy_chromosome_cpgs = nrow(off_policy),
    n_vmrs_with_zero_tested = denom[n_tested_cpgs == 0, .N],
    n_donors = length(donors_seen)
)
write_atomic(summary_dt, file.path(out, "cpg-set-summary.tsv"))
write_atomic(excl[, .N, by = reason], file.path(out, "exclusion-reasons.tsv"))
print(summary_dt)
print(excl[, .N, by = reason])

## Member CpGs must equal tested plus excluded plus those absent from every
## matrix. A remainder means CpGs vanished silently.
unaccounted <- nrow(membership) - length(tested) - nrow(excl)
if (unaccounted != 0L) {
    stop("CpG accounting does not balance: ", nrow(membership), " member, ",
         length(tested), " tested, ", nrow(excl), " excluded, ",
         unaccounted, " unaccounted.")
}

append_manifest(list(dir = run_dir), list(
    n_donors = length(donors_seen),
    donor_checksum = donor_checksum(donors_seen),
    n_member_cpgs = nrow(membership),
    n_tested_cpgs = length(tested)
))
if (!is.null(smoke_chroms)) {
    append_manifest(list(dir = run_dir), list(smoke_chroms = smoke_chroms))
}
