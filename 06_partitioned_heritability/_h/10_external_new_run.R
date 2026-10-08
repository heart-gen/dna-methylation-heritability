#!/usr/bin/env Rscript
#### 06_partitioned_heritability -- open an external-annotation benchmark run ####
##
## Usage:
##   Rscript _h/10_external_new_run.R --cohort AA [--allow-unlocked] [--run-id ID]
##
## Opens `sldsc-{cohort}-external-{date}` and stages everything the LD-score and
## regression stages read, so the run is self-contained:
##
##   annotation/<NAME>.hg19.bed      one binary annotation per row of
##                                   config/sldsc_external_annotations.yml with
##                                   status ready, plus VMR_TESTED_<REGION>
##   annotation/annotation-list.tsv  what was staged, and what was skipped and why
##   sumstats/<trait>.sumstats.gz    the munged files of the accepted VMR runs
##
## "Exactly the same pipeline" is enforced, not described:
##   - the three accepted runs named in config must be the accepted Module 06 runs;
##   - their munged sumstats must have identical CONTENT across cells (the gzip
##     bytes differ by timestamp, so the check is on the decompressed stream),
##     and those are the files copied here;
##   - the LD reference block is read from config/partitioned_heritability.yml at
##     its ld_reference_arm by the later stages, and its checksum is recorded.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
})

MODULE <- "06_partitioned_heritability"
MODULE_TAG <- "sldsc"

opts <- parse_v2_args(require = "cohort")
allow_unlocked <- isTRUE(opts$allow_unlocked)
cohort <- opts$cohort

ph <- load_config("partitioned_heritability")
ext <- load_config("sldsc_external_annotations")
assert_locked(list(partitioned_heritability = ph, sldsc_external_annotations = ext),
              allow_unlocked = allow_unlocked)
arm <- ph$ld_reference_arm
if (is.null(ph$ld_references[[arm]])) stop("ld_reference_arm '", arm, "' has no entry")

regions <- names(ext$accepted_vmr_runs)
up <- lapply(regions, function(re)
    require_accepted_upstream(MODULE, cohort, re, allow_unaccepted = allow_unlocked))
names(up) <- regions
for (re in regions) {
    if (!identical(as.character(up[[re]]$run_id), ext$accepted_vmr_runs[[re]])) {
        stop("config names ", ext$accepted_vmr_runs[[re]], " for ", re,
             " but the accepted Module 06 run is ", up[[re]]$run_id)
    }
}
vmr_dir <- function(re) file.path(repo_root(), MODULE, "_m", "runs", up[[re]]$run_id)

## --------------------------------------------------- identical GWAS inputs
traits <- vapply(ph$traits, `[[`, character(1), "name")
content_md5 <- function(f) {
    out <- system2("bash", c("-c", shQuote(paste("zcat", shQuote(f), "| md5sum"))),
                   stdout = TRUE)
    sub(" .*", "", out)
}
sumstat_check <- rbindlist(lapply(traits, function(t) {
    md5 <- vapply(regions, function(re)
        content_md5(file.path(vmr_dir(re), "sumstats", paste0(t, ".sumstats.gz"))),
        character(1))
    data.table(trait = t, region = regions, content_md5 = md5,
               identical_across_cells = length(unique(md5)) == 1L)
}))
if (!all(sumstat_check$identical_across_cells)) {
    print(sumstat_check[identical_across_cells == FALSE])
    stop("The accepted runs' munged sumstats differ in content; there is no single ",
         "GWAS set to reuse.")
}

## ------------------------------------------------------------------ the run
if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked)")
}
run <- new_run(
    module = MODULE_TAG, cohort = cohort, region = "external",
    module_root = file.path(repo_root(), MODULE), run_id = opts$run_id,
    upstream = stats::setNames(lapply(regions, function(re) up[[re]]$run_id),
                               paste0("partitioned_heritability_", regions)),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        analysis = "external_annotation_benchmark",
        gating = "FALSE",
        config_partitioned_heritability_sha256 = attr(ph, "config_sha256"),
        config_sldsc_external_annotations_sha256 = attr(ext, "config_sha256"),
        ld_reference_arm = arm,
        ld_reference_label = ph$ld_references[[arm]]$label,
        n_traits = length(traits),
        sumstats_source_run = up[[regions[1]]]$run_id
    )
)
for (d in c("annotation", "sumstats", "ldscores", "results", "code/config")) {
    dir.create(file.path(run$dir, d), recursive = TRUE, showWarnings = FALSE)
}
invisible(file.copy(file.path(repo_root(), "config",
                              c("partitioned_heritability.yml", "sldsc_external_annotations.yml")),
                    file.path(run$dir, "code", "config")))
for (t in traits) {
    ok <- file.copy(file.path(vmr_dir(regions[1]), "sumstats", paste0(t, ".sumstats.gz")),
                    file.path(run$dir, "sumstats", paste0(t, ".sumstats.gz")))
    if (!ok) stop("could not stage sumstats for ", t)
}
write_atomic(sumstat_check, file.path(run$dir, "results", "sumstats-identity-check.tsv"))

## ------------------------------------------------------------ annotations
AUTOSOMES <- paste0("chr", 1:22)
write_bed <- function(gr, name) {
    gr <- GenomicRanges::reduce(gr[as.character(seqnames(gr)) %in% AUTOSOMES],
                                ignore.strand = TRUE)
    gr <- sort(gr)
    f <- file.path(run$dir, "annotation", paste0(name, ".hg19.bed"))
    fwrite(data.table(as.character(seqnames(gr)), start(gr) - 1L, end(gr)), f,
           sep = "\t", col.names = FALSE)
    list(n = length(gr), bp = sum(as.numeric(width(gr))), sha = file_sha256(f))
}
read_bed3 <- function(f) {
    dt <- fread(cmd = paste(if (grepl("\\.gz$", f)) "zcat" else "cat", shQuote(f)),
                header = FALSE, select = 1:3, col.names = c("chrom", "start", "end"))
    GRanges(dt$chrom, IRanges(dt$start + 1L, dt$end))
}

listing <- list()
for (a in ext$annotations) {
    status <- a$status %||% "ready"
    if (!identical(status, "ready") || is.null(a$bed)) {
        listing[[length(listing) + 1]] <- data.table(
            annotation = a$name, role = a$role, label = a$label, status = status,
            source = a$source_url %||% NA_character_)
        next
    }
    src <- file.path(repo_root(), a$bed)
    if (!file.exists(src)) stop("Missing ", src, "; run inputs/supportfiles/_h/03_build_rizzardi_dmr_asset.py")
    w <- write_bed(read_bed3(src), a$name)
    listing[[length(listing) + 1]] <- data.table(
        annotation = a$name, role = a$role, label = a$label, status = "ready",
        source = a$source_url, n_intervals = w$n, bp = w$bp, bed_sha256 = w$sha)
}
if (isTRUE(ext$vmr_membership_comparators)) {
    for (re in regions) {
        nm <- paste0("VMR_TESTED_", toupper(re))
        w <- write_bed(read_bed3(file.path(vmr_dir(re), "annotation", "annotation-hg19.bed")), nm)
        listing[[length(listing) + 1]] <- data.table(
            annotation = nm, role = "vmr_membership_standalone",
            label = paste("Tested VMR universe,", re, "(accepted", up[[re]]$run_id, ")"),
            status = "ready", source = up[[re]]$run_id,
            n_intervals = w$n, bp = w$bp, bed_sha256 = w$sha)
    }
}
listing <- rbindlist(listing, fill = TRUE)
## Names must not nest: the regression stage identifies its .results row by name.
nm <- listing[status == "ready", annotation]
for (x in nm) for (y in nm) if (x != y && startsWith(y, x)) stop("Annotation names nest: ", x, " / ", y)
write_atomic(listing, file.path(run$dir, "annotation", "annotation-list.tsv"))
print(listing)
message("[06] external run ", run$run_id, " opened")
cat(run$run_id, "\n", sep = "")
