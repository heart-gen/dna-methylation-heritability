#### 10 / QC: cross-region donor concordance and sample-swap detection ####
##
## PI decision D3 (2026-08-26). Migrates the one QC check in
## qc_analysis/_h/04.qc_tissue_compare.qmd that Module 01 does not perform:
## Module 01 does coverage and methylation QC within each cohort x region cell
## independently, and so can never notice that a donor's caudate and DLPFC
## samples disagree. In a three-region post-mortem design a mislabelled tube is
## a live risk, and it is invisible to every per-cell check in the pipeline.
##
## The test. For donors sampled in more than one region, methylation profiles
## at a common set of CpGs should correlate more strongly with the SAME donor's
## other region than with a typical OTHER donor. The reported flag is
## z_self < 0: the donor resembles the average other donor more than itself.
##
## Power caveat, stated up front. Module 01 masks C->T SNP CpGs
## (00_prepare.R remove_ct_snps), which are exactly the genotype-driven sites
## a methylation fingerprint would rely on. What remains is a real but weak
## identity signal -- 45-65% of donors are their own single best match out of
## ~150 candidates, against a ~0.7% chance rate. That is far above chance but
## not clean enough to adjudicate an individual donor, so this is a SCREEN
## that nominates candidates for follow-up, not a verdict. A definitive check
## needs the unmasked C->T CpGs.
##
## This is deliberately not a genotype concordance check: the genotypes come
## from one blood/tissue draw per donor and are shared across the donor's
## regions by construction, so they cannot detect a swap introduced at the
## methylation assay. Methylation is where the swap would show.
##
## Scope. One autosome (default chr22, the smallest) is enough: the statistic
## is a donor-level correlation over thousands of CpGs, not a locus result.
## Nothing downstream consumes this; it is a supplementary QC figure.
##
## Outputs
##   tables/qc_cross_region_concordance.tsv   per donor x region pair
##   tables/qc_swap_candidates.tsv            donors failing the check
##   figures/figureS_sample_integrity         same- vs different-donor r
##
## Usage:
##   Rscript 05_qc_sample_integrity.R --run-id ID [--cohort AA] [--chrom 22]
##                                    [--max-cpgs 5000]

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "10_integrated_manuscript_outputs", "_h",
                 "00_figure_theme.R"))

opts    <- parse_v2_args(require = c("run_id"))
cohorts <- load_config("cohorts")
regions <- cohorts$regions
cohort  <- if (is.null(opts$cohort)) cohorts$primary else opts$cohort
chrom   <- if (is.null(opts$chrom)) "22" else as.character(opts$chrom)
max_cpg <- if (is.null(opts$max_cpgs)) 5000L else as.integer(opts$max_cpgs)

module_root <- file.path(V2_ROOT, "10_integrated_manuscript_outputs")
run_dir   <- file.path(module_root, "_m", "runs", opts$run_id)
fig_dir   <- file.path(run_dir, "figures")
data_dir  <- file.path(run_dir, "source_data")
table_dir <- file.path(run_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

SCRIPT <- "10_integrated_manuscript_outputs/_h/05_qc_sample_integrity.R"

## --------------------------------------------------------------- locate data
runs <- vapply(regions, function(r) {
    require_accepted_upstream("01_vmr_catalog", cohort = cohort, region = r)$run_id
}, "")
names(runs) <- regions

meth_file <- function(r) file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs",
                                   runs[[r]], "cpg", paste0("chr_", chrom),
                                   "cpg_meth.phen")
for (r in regions) if (!file.exists(meth_file(r)))
    stop("no chr", chrom, " methylation matrix for ", r, ": ", meth_file(r))

## ------------------------------------------------- common CpGs across regions
##
## The matrices are donors x CpGs with position names in the header and are
## ~1.3 GB each, so the header alone is read first and only the selected
## columns are then streamed through `cut`. Reading whole matrices to keep
## 5,000 columns would cost ~4 GB of I/O for no gain.

header_of <- function(r) {
    h <- strsplit(readLines(meth_file(r), n = 1L), "\t", fixed = TRUE)[[1]]
    h[-(1:2)]                                   # drop FID, IID
}
pos <- lapply(regions, header_of); names(pos) <- regions
common <- Reduce(intersect, pos)
if (length(common) < 100)
    stop("only ", length(common), " CpGs shared across regions on chr", chrom)

## Evenly spaced rather than random: spreads the pool across the chromosome and
## makes the selection reproducible without carrying a seed. A pool larger than
## the final set is read because the informative CpGs are chosen by variance
## below.
pool_n <- min(max(max_cpg * 4L, 20000L), length(common))
pool <- common[unique(round(seq(1, length(common), length.out = pool_n)))]
message("[qc] ", length(common), " common chr", chrom, " CpGs; pool ", length(pool))

## Column selection goes through awk with the indices in a file, not `cut -f`
## with them on the command line: a pool of 20,000 indices overruns the
## argument limit. The matrices are only a few hundred rows, so this is cheap.
read_region <- function(r) {
    idx <- match(pool, pos[[r]]) + 2L           # +2 for FID, IID
    idx_file <- tempfile(fileext = ".idx")
    on.exit(unlink(idx_file), add = TRUE)
    writeLines(as.character(idx), idx_file)
    prog <- paste0(
        'BEGIN{while((getline l < "', idx_file, '")>0) k[++n]=l+0}',
        '{s=$1; for(i=1;i<=n;i++) s=s"\t"$(k[i]); print s}')
    d <- fread(cmd = paste("awk", shQuote(prog), shQuote(meth_file(r))),
               header = TRUE, colClasses = list(character = 1))
    setnames(d, 1L, "brnum")
    m <- as.matrix(d[, -1L]); rownames(m) <- d$brnum
    colnames(m) <- pool
    m
}
mats <- lapply(regions, read_region); names(mats) <- regions

## Why centering is not optional. Raw cross-region correlation is dominated by
## the CpG-level mean profile, which every donor shares: uncentered, r(self)
## and r(other) both sit around 0.97 and the check cannot separate them. What
## identifies a donor is the DEVIATION from the CpG mean, so each CpG is
## centred within its region before any correlation is taken.
##
## Then keep the most variable CpGs: a CpG that barely varies between donors
## carries no identity information and only adds noise to the correlation.
mats <- lapply(mats, function(m) sweep(m, 2L, colMeans(m, na.rm = TRUE), "-"))

cpg_var <- Reduce(`+`, lapply(mats, function(m)
    matrixStats::colVars(m, na.rm = TRUE))) / length(mats)
keep <- order(cpg_var, decreasing = TRUE)[seq_len(min(max_cpg, length(pool)))]
sel  <- pool[sort(keep)]
mats <- lapply(mats, function(m) m[, sort(keep), drop = FALSE])
message("[qc] using the ", length(sel), " most variable of them")

## ------------------------------------------------------- pairwise concordance
##
## For each ordered region pair and each donor present in both, r(same donor)
## against the best r to any OTHER donor. Correlations use complete pairs only.

pairs <- combn(regions, 2, simplify = FALSE)

concordance <- rbindlist(lapply(pairs, function(pr) {
    a <- mats[[pr[1]]]; b <- mats[[pr[2]]]
    shared <- intersect(rownames(a), rownames(b))
    if (length(shared) < 3) return(NULL)
    ## cor() on the transposed matrices gives donor x donor across CpGs.
    cm <- suppressWarnings(stats::cor(t(a[shared, , drop = FALSE]),
                                      t(b[shared, , drop = FALSE]),
                                      use = "pairwise.complete.obs"))
    self  <- diag(cm)
    other <- cm; diag(other) <- NA_real_
    best_other <- apply(other, 1, max, na.rm = TRUE)
    best_match <- shared[apply(other, 1, function(x) which.max(replace(x, is.na(x), -Inf)))]
    mu <- rowMeans(other, na.rm = TRUE)
    sdv <- apply(other, 1, stats::sd, na.rm = TRUE)
    ## Rank of the donor's own other-region sample among all candidates,
    ## 1 = the donor is its own best match.
    rk <- vapply(seq_along(shared), function(i)
        1L + sum(other[i, ] > self[i], na.rm = TRUE), integer(1))
    data.table(region_a = pr[1], region_b = pr[2], brnum = shared,
               n_candidates = length(shared),
               r_self = as.numeric(self),
               r_mean_other = as.numeric(mu),
               r_best_other = as.numeric(best_other),
               best_other_donor = best_match,
               rank_self = rk,
               z_self = as.numeric((self - mu) / sdv))
}))

if (nrow(concordance) == 0)
    stop("no donors are shared between any two regions; the swap check cannot run")

## Calibration. "Self must be the single best match" is the wrong criterion on
## this matrix: Module 01 masks C->T SNP CpGs (00_prepare.R remove_ct_snps),
## which are exactly the genotype-driven sites that carry donor identity, so
## the residual fingerprint is real but weak. Empirically 45-65% of donors are
## their own top match out of ~150 candidates -- vastly above the ~0.7% chance
## rate, so there is signal -- but requiring rank 1 would flag a third of a
## healthy cohort.
##
## The reported flag is therefore the directional one: a donor whose own
## other-region sample resembles it LESS than a typical other donor does
## (z_self < 0). That is the pattern an actual swap produces.
concordance[, swap_candidate := z_self < 0]
setorder(concordance, z_self)

write_atomic(concordance, file.path(table_dir, "qc_cross_region_concordance.tsv"))
write_atomic(concordance[swap_candidate == TRUE],
             file.path(table_dir, "qc_swap_candidates.tsv"))

n_bad <- concordance[, sum(swap_candidate)]
message("[qc] ", nrow(concordance), " donor x region-pair comparisons; ",
        n_bad, " swap candidate(s)")
if (n_bad > 0)
    warning("SWAP CANDIDATES in ", cohort, ": see qc_swap_candidates.tsv. ",
            "These donors resemble a typical other donor across regions MORE ",
            "than they resemble themselves (z_self < 0). Resolve before ",
            "submission. Note this screen is underpowered by design: the ",
            "C->T SNP CpGs that carry donor identity are masked upstream.")

## ------------------------------------------------------------------- figure
plot_dt <- melt(concordance,
                id.vars = c("region_a", "region_b", "brnum"),
                measure.vars = c("r_self", "r_mean_other", "r_best_other"),
                variable.name = "comparison", value.name = "r")
plot_dt[, comparison := factor(
    comparison, levels = c("r_self", "r_mean_other", "r_best_other"),
    labels = c("Same donor", "Mean other donor", "Best other donor"))]
plot_dt[, pair := paste0(REGION_LABELS[region_a], " vs ", REGION_LABELS[region_b])]

p <- ggplot(plot_dt, aes(comparison, r, colour = comparison)) +
    geom_line(aes(group = brnum), colour = PAL_NULL, alpha = 0.35,
              linewidth = 0.25) +
    geom_point(size = 0.7, alpha = 0.75) +
    facet_wrap(~ pair, nrow = 1) +
    scale_colour_manual(values = c("Same donor" = PAL_BLUE,
                                   "Mean other donor" = PAL_NULL,
                                   "Best other donor" = PAL_RUST),
                        guide = "none") +
    labs(x = NULL,
         y = paste0("Cross-region correlation (chr", chrom, ", ",
                    length(sel), " CpGs)")) +
    BASE_THEME + NO_TITLES +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_figure(p, "figureS_sample_integrity", FIG_WIDTH_FULL, 3.2, fig_dir)
write_source_data(
    concordance, "figureS_sample_integrity",
    source_run_id = unname(runs),
    source_table  = paste0("01_vmr_catalog cpg/chr_", chrom, "/cpg_meth.phen"),
    script        = SCRIPT,
    filter_desc   = paste0(cohort, " donors present in both regions of each pair; ",
                           length(sel), " evenly spaced CpGs common to all regions"),
    data_dir      = data_dir)

message("[10] sample-integrity QC written to ", run_dir)
