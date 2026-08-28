#### inputs/supportfiles / 01_build_array_universe: array probe universes ####
##
## Builds the assayed-CpG universe for an Illumina methylation array, lifted to
## hg38, in the schema the 450K universe already uses:
##
##     probe_id    chrom    pos_1based
##
## Why this exists: AGENTS.md 2.2 claims "WGBS coverage of methylation outside
## array-accessible CpGs". That is a comparison against a platform we did not
## run, so it needs the array's probe coordinates as a reference annotation --
## the same category as repeat-masker or the ENCODE blacklist. No array
## intensities, samples, or methylation values are involved; only coordinates.
##
## The 450K universe was first built inline by
## meqtl-validation/03_external_meqtl_validation/_h/03_harmonize_external_meqtls.py
## (ensure_450k_universe). This script generalizes that logic so Module 01 can
## compare against both 450K and EPIC without depending on a meQTL script.
##
## Usage:
##   Rscript 01_build_array_universe.R --platform EPIC
##   Rscript 01_build_array_universe.R --platform 450K

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

opts <- parse_v2_args(require = "platform")
platform <- toupper(opts$platform)

## Only probes that lift cleanly are kept: an unlifted probe has no hg38
## coordinate, so it can neither match nor exclude a WGBS CpG.
PLATFORMS <- list(
    "450K" = list(
        pkg = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        out = "450k_universe_hg38.tsv.gz",
        expected = 485512L),
    "EPIC" = list(
        pkg = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
        out = "epic_universe_hg38.tsv.gz",
        expected = 865859L)
)

if (!platform %in% names(PLATFORMS)) {
    stop("Unknown --platform ", platform, "; expected one of ",
         paste(names(PLATFORMS), collapse = ", "))
}
spec <- PLATFORMS[[platform]]

support_dir <- file.path(V2_ROOT, "inputs", "supportfiles", "_m")
out_path <- file.path(support_dir, spec$out)

CHAIN <- file.path(V2_ROOT, "meqtl-validation", "03_external_meqtl_validation",
                   "_m", "support", "hg19ToHg38.over.chain.gz")
LIFTOVER <- "/projects/p32505/opt/envs/genomics/bin/liftOver"

for (dep in c(CHAIN, LIFTOVER)) {
    if (!file.exists(dep)) stop("Required dependency not found: ", dep)
}
if (!requireNamespace(spec$pkg, quietly = TRUE)) {
    stop("Annotation package not installed: ", spec$pkg,
         "\nInstall with BiocManager::install() in the epigenomics env.")
}

## -------------------------------------------------------- manifest -> hg19 BED
suppressPackageStartupMessages(library(spec$pkg, character.only = TRUE))
anno <- minfi::getAnnotation(get(spec$pkg))

hg19 <- data.table(
    chr   = as.character(anno$chr),
    start = as.integer(anno$pos) - 1L,
    end   = as.integer(anno$pos),
    name  = rownames(anno))
hg19 <- hg19[!is.na(chr) & !is.na(start)]
n_source <- nrow(hg19)
message("[", platform, "] manifest probes: ", n_source)

if (n_source < 0.95 * spec$expected) {
    stop("Manifest has ", n_source, " probes, far below the expected ~",
         spec$expected, " for ", platform, "; wrong annotation package?")
}

tmp <- tempfile(pattern = paste0(platform, "_universe_"))
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
bed_hg19 <- file.path(tmp, "hg19.bed")
bed_hg38 <- file.path(tmp, "hg38.bed")
bed_unmapped <- file.path(tmp, "unmapped.bed")

fwrite(hg19, bed_hg19, sep = "\t", col.names = FALSE, quote = FALSE)

## ------------------------------------------------------------------ liftOver
status <- system2(LIFTOVER, c(bed_hg19, CHAIN, bed_hg38, bed_unmapped))
if (status != 0L) stop("liftOver failed with exit status ", status)

lifted <- fread(bed_hg38, header = FALSE,
                col.names = c("chrom", "start0", "pos_1based", "probe_id"),
                colClasses = list(character = c(1, 4)))
lifted <- unique(lifted, by = "probe_id")

frac <- nrow(lifted) / n_source
message("[", platform, "] lifted ", nrow(lifted), "/", n_source,
        " probes to hg38 (", sprintf("%.2f%%", 100 * frac), ")")
if (frac < 0.95) {
    stop("Only ", sprintf("%.2f%%", 100 * frac), " of probes lifted; ",
         "refusing to write a truncated universe.")
}

universe <- lifted[, .(probe_id, chrom, pos_1based)]
setkey(universe, chrom, pos_1based)

fwrite(universe, out_path, sep = "\t", compress = "gzip")
message("[done] wrote ", nrow(universe), " probes to ", out_path)

## ------------------------------------------------- asset manifest registration
##
## Printed rather than appended: annotation_asset_manifest.tsv is edited once,
## deliberately, with the checksum verified by eye.
sha256 <- system2("sha256sum", shQuote(out_path), stdout = TRUE)
sha256 <- sub(" .*$", "", sha256)
message("\nAdd this row to inputs/supportfiles/_m/annotation_asset_manifest.tsv:\n")
cat(paste(spec$out, out_path,
          sprintf("Illumina %s manifest lifted to hg38 (%d/%d probes)",
                  platform, nrow(universe), n_source),
          file.size(out_path), sha256,
          format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), "ready",
          sprintf("array probe universe used by 01_vmr_catalog/04_turnover.R to quantify WGBS coverage outside %s-accessible CpGs (AGENTS.md 2.2)", platform),
          sep = "\t"), "\n")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
