#### 10 / Table 1: donor demographics, and the ancestry-confirmation panel ####
##
## PI decision D3 (2026-08-26): Module 10 owns Table 1 and the cohort QC
## figures, reporting the AA primary arm and the all_individuals sensitivity
## arm side by side.
##
## Why this exists at all. The legacy Table 1
## (sample_summary/_h/01.fancy_table.sample_summary.R) derived its donor set
## from vmr-analysis/all_individuals/{region}/_m/samples.txt, which defect V1
## invalidated, and it honoured the legacy sample blacklists that v2 retired.
## Eight donors (Br1249 Br1303 Br1371 Br1552 Br1693 Br1700 Br1883 Br1927) are
## therefore missing from the published table but present in v2. The table is
## not merely unmigrated; it reports the wrong cohort. See config/cohorts.yml
## and MIGRATION_MANIFEST.tsv.
##
## The donor set here is read from the ACCEPTED Module 01 catalog runs
## (vmr/donors_plink.txt), never from a phenotype file, so the table can only
## ever describe donors that actually entered the analysis.
##
## gtsummary and gt are not in the epigenomics environment, so the formatted
## table is assembled with data.table and written as TSV (source data) plus a
## LaTeX booktabs fragment for the manuscript.
##
## Outputs
##   tables/table1_cohort.tsv        long-form, one row per statistic
##   tables/table1_cohort.tex        booktabs fragment, both arms
##   figures/figureS_ancestry_pcs    snpPC1/2, donors over 1000 Genomes
##
## Usage:
##   Rscript 04_table1_cohort.R --run-id fig-all-20260827

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "10_integrated_manuscript_outputs", "_h",
                 "00_figure_theme.R"))

suppressPackageStartupMessages(library(scales))

opts    <- parse_v2_args(require = c("run_id"))
cohorts <- load_config("cohorts")
regions <- cohorts$regions
arms    <- cohorts$arms
paths   <- load_config("paths")

module_root <- file.path(V2_ROOT, "10_integrated_manuscript_outputs")
run_dir   <- file.path(module_root, "_m", "runs", opts$run_id)
fig_dir   <- file.path(run_dir, "figures")
data_dir  <- file.path(run_dir, "source_data")
table_dir <- file.path(run_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

SCRIPT <- "10_integrated_manuscript_outputs/_h/04_table1_cohort.R"

## --------------------------------------------------------- accepted donors
##
## require_accepted_upstream() is the gate (AGENTS.md 6); reading the run
## directory without it would let this table describe an unaccepted catalog.
## The gate returns the accepted row, so the run ID comes from the README of
## record rather than being reconstructed from a naming convention here.

accepted <- new.env(parent = emptyenv())
CATALOG_RUN <- function(arm, r) get(paste(arm, r), envir = accepted)

donors <- rbindlist(lapply(arms, function(arm) {
    rbindlist(lapply(regions, function(r) {
        up <- require_accepted_upstream("01_vmr_catalog", cohort = arm, region = r)
        assign(paste(arm, r), up$run_id, envir = accepted)
        f <- file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs",
                       up$run_id, "vmr", "donors_plink.txt")
        d <- fread(f, header = FALSE, col.names = c("brnum", "iid"),
                   colClasses = "character")
        d[, `:=`(arm = arm, region = r)][]
    }))
}))

## Guard: the accepted donor count must equal the locked design_n. A silent
## mismatch here is exactly the failure mode that produced the wrong legacy
## table, so it is fatal rather than a warning.
for (a in arms) for (r in regions) {
    want <- cohorts$donor_counts[[a]][[r]]$design_n
    got  <- donors[arm == a & region == r, .N]
    if (!is.null(want) && !is.na(want) && want != got)
        stop("design_n mismatch for ", a, "/", r, ": locked ", want,
             ", accepted run has ", got)
}

## ------------------------------------------------------------- phenotypes
## NOT paths$phenotype_table -- that key still points at the legacy
## sample_summary/_m/phenotype_data.tsv this script exists to replace.
## cohort_def() resolves the one table both arms share (config/cohorts.yml).
pheno_file <- cohort_def(arms[1])$phenotype_table
for (a in arms) if (!identical(cohort_def(a)$phenotype_table, pheno_file))
    stop("arms disagree on phenotype_table; this table assumes one shared source")
pheno <- fread(pheno_file)
pheno[, brnum := as.character(brnum)]

## One phenotype row per donor x region. Join on both so a donor sampled in
## several regions contributes its own row to each region's column.
dt <- merge(donors, pheno, by = c("brnum", "region"), all.x = TRUE)
if (anyNA(dt$agedeath))
    stop(sum(is.na(dt$agedeath)), " accepted donors have no phenotype row; ",
         "the table would silently under-report them")

dt[, `:=`(
    sex       = factor(sex, levels = c("F", "M"), labels = c("Female", "Male")),
    primarydx = factor(as.character(primarydx)),
    race      = factor(race, levels = c("AA", "EA"),
                       labels = c("Black American", "Non-Hispanic white American")),
    region_lab = as_region(region),
    arm_lab    = factor(unname(COHORT_LABELS[arm]),
                        levels = unname(COHORT_LABELS[arms])))]

## ------------------------------------------------------------ Table 1 body
##
## Long form (one row per statistic) rather than a wide pre-rendered grid: the
## source-data convention wants the numbers machine-readable, and the LaTeX
## fragment is generated from this, so the table and its source cannot drift.

fmt_n_pct <- function(n, d) sprintf("%d (%.1f%%)", n, 100 * n / d)

summarise_cell <- function(d) {
    n <- nrow(d)
    out <- list(
        data.table(variable = "N", level = "", value = as.character(n)),
        data.table(variable = "Age at death, years",
                   level = "mean (SD)",
                   value = sprintf("%.1f (%.1f)", mean(d$agedeath),
                                   stats::sd(d$agedeath))),
        data.table(variable = "Age at death, years",
                   level = "range",
                   value = sprintf("%.1f-%.1f", min(d$agedeath), max(d$agedeath)))
    )
    for (v in c("sex", "primarydx", "race")) {
        lv <- d[[v]]
        if (all(is.na(lv))) next
        tab <- table(droplevels(lv))
        if (length(tab) < 2 && v == "race") next   # AA arm is single-race
        lab <- c(sex = "Sex", primarydx = "Primary diagnosis",
                 race = "Donor group")[[v]]
        out <- c(out, list(data.table(
            variable = lab, level = names(tab),
            value = vapply(seq_along(tab),
                           function(i) fmt_n_pct(tab[[i]], nrow(d)), ""))))
    }
    if (!all(is.na(d$pmi)))
        out <- c(out, list(data.table(
            variable = "Post-mortem interval, hours", level = "mean (SD)",
            value = sprintf("%.1f (%.1f)", mean(d$pmi, na.rm = TRUE),
                            stats::sd(d$pmi, na.rm = TRUE)))))
    rbindlist(out)
}

table1 <- dt[, summarise_cell(.SD), by = .(arm, arm_lab, region, region_lab)]
setcolorder(table1, c("arm", "arm_lab", "region", "region_lab",
                      "variable", "level", "value"))

write_source_data(
    table1, "table1_cohort",
    source_run_id = unlist(lapply(arms, function(a) vapply(regions, function(r) CATALOG_RUN(a, r), ""))),
    source_table  = "01_vmr_catalog vmr/donors_plink.txt + inputs/phenotypes/_m/phenotypes-all.tsv",
    script        = SCRIPT,
    filter_desc   = "donors in the accepted Module 01 catalog run for each arm x region",
    data_dir      = data_dir)
write_atomic(table1, file.path(table_dir, "table1_cohort.tsv"))

## ------------------------------------------------------------ LaTeX output
esc <- function(x) gsub("%", "\\\\%", gsub("&", "\\\\&", x))

tex_for_arm <- function(a) {
    d <- table1[arm == a]
    rows <- unique(d[, .(variable, level)])
    hdr <- paste(c("Characteristic", "",
                   as.character(unique(d$region_lab))), collapse = " & ")
    body <- vapply(seq_len(nrow(rows)), function(i) {
        v <- rows$variable[i]; l <- rows$level[i]
        vals <- vapply(regions, function(r) {
            x <- d[variable == v & level == l & region == r, value]
            if (length(x) == 0) "--" else x[1]
        }, "")
        paste(c(esc(v), esc(l), esc(vals)), collapse = " & ")
    }, "")
    c(sprintf("%% %s", unname(COHORT_LABELS[a])),
      "\\begin{tabular}{ll" , paste(rep("r", length(regions)), collapse = ""), "}",
      "\\toprule", paste0(hdr, " \\\\"), "\\midrule",
      paste0(body, " \\\\"), "\\bottomrule", "\\end{tabular}", "")
}

writeLines(unlist(lapply(arms, tex_for_arm)),
           file.path(table_dir, "table1_cohort.tex"))
message("[table] table1_cohort.tsv + .tex (", nrow(table1), " statistic rows)")

## ------------------------------------------- Figure S: ancestry confirmation
##
## Confirms the recorded donor group against genotype, with 1000 Genomes as the
## external frame of reference. This is a QC panel, not an ancestry-biology
## claim: AGENTS.md 2.3 forbids attributing differences to ancestry-specific
## biology, and nothing downstream consumes this figure.

ref_pc  <- fread(paths$reference_1kgp$eigenvec)
setnames(ref_pc, c("IID", paste0("snpPC", 1:10))[seq_len(ncol(ref_pc))])
## The panel file carries trailing empty header fields; name only what fread
## actually returns rather than letting it fill.
ref_pop <- fread(paths$reference_1kgp$sample_panel, fill = TRUE)
setnames(ref_pop, seq_len(4L), c("IID", "pop", "super_pop", "gender"))
ref <- merge(ref_pc[, .(IID, snpPC1, snpPC2)], ref_pop[, .(IID, super_pop)],
             by = "IID")

## One point per donor (not per donor x region), from the widest arm.
wide_arm <- if ("all_individuals" %in% arms) "all_individuals" else arms[1]
donor_pc <- unique(dt[arm == wide_arm, .(brnum, race, snpPC1, snpPC2)])
donor_pc <- donor_pc[!is.na(snpPC1)]

pcs <- ggplot() +
    geom_point(data = ref, aes(snpPC1, snpPC2, colour = super_pop),
               size = 0.6, alpha = 0.45) +
    geom_point(data = donor_pc, aes(snpPC1, snpPC2, shape = race),
               colour = PAL_CHARCOAL, size = 1.5, stroke = 0.5, fill = NA) +
    scale_colour_brewer(palette = "Set2", name = "1000 Genomes") +
    scale_shape_manual(values = c(21, 24), name = "This study") +
    labs(x = "Genotype PC1", y = "Genotype PC2") +
    BASE_THEME + NO_TITLES

save_figure(pcs, "figureS_ancestry_pcs", FIG_WIDTH_THREEQ, 3.4, fig_dir)
write_source_data(
    rbind(ref[, .(source = "1000 Genomes", group = super_pop, snpPC1, snpPC2)],
          donor_pc[, .(source = "this study", group = as.character(race),
                       snpPC1, snpPC2)]),
    "figureS_ancestry_pcs",
    source_run_id = vapply(regions, function(r) CATALOG_RUN(wide_arm, r), ""),
    source_table  = "phenotypes-all.tsv snpPC1-2 + 1kGP-pc.eigenvec",
    script        = SCRIPT,
    filter_desc   = paste0("unique donors in the accepted ", wide_arm,
                           " catalog runs with non-missing genotype PCs"),
    data_dir      = data_dir)

message("[10] Table 1 and ancestry panel written to ", run_dir)
