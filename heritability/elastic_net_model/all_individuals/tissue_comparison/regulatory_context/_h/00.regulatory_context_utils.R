#### Shared utilities for VMR regulatory-context analyses ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(GenomicRanges)
  library(SummarizedExperiment)
})

TISSUES <- c("caudate", "dlpfc", "hippocampus")
H2_GROUP_LEVELS <- c("Heritable", "Non-heritable", "Low prediction")

## Matched shared VMRs use one genomic list (matched_r2_0.3.tsv) with per-population
## h2_unscaled / r_squared_cv columns and a single concordant h2_category column.
## Association tests and genomic proximity do not depend on population; only
## continuous-h2 quintile summaries differ by population.
SHARED_VMR_CANONICAL_POPULATION <- "AA"

rna_file_map <- list(
  expression = c(
    caudate = "rse-gene.bsp3.caudate-n487.gencode-v47.RData",
    dlpfc = "rse-gene.bsp2.dlpfc-n500.gencode-v47.RData",
    hippocampus = "rse-gene.bsp2.hippocampus-n452.gencode-v47.RData"
  ),
  psi = c(
    caudate = "rse-psi.bsp3.caudate-n487.gencode-v47.RData",
    dlpfc = "rse-psi.bsp2.dlpfc-n500.gencode-v47.RData",
    hippocampus = "rse-psi.bsp2.hippocampus-n452.gencode-v47.RData"
  )
)

message2 <- function(...) {
  message(sprintf(...))
}

## Output subfolders under regulatory_context/_m/{tissue}/{pop}/expression/ from
## 01.run_local_associations.R (ABC vs nearest-gene passes).
nearest_gene_run_tag <- function(window_bp = 250000L) {
  paste0("nearest_gene_window_", as.integer(window_bp) / 1000, "kb")
}

resolve_regulatory_run_tags <- function(modality, window_bp = 250000L,
                                        expression_layers = "both") {
  if (modality == "psi") {
    return(paste0("window_", as.integer(window_bp) / 1000, "kb"))
  }
  if (modality != "expression") {
    stop("resolve_regulatory_run_tags supports modality expression or psi")
  }
  el <- trimws(tolower(as.character(expression_layers)))
  ng <- nearest_gene_run_tag(window_bp)
  if (el %in% c("both", "all", "")) {
    return(c("abc", ng))
  }
  if (el == "abc") {
    return("abc")
  }
  if (el %in% c("nearest_gene", "nearest", "ng")) {
    return(ng)
  }
  stop(
    "expression_link_layers must be both, abc, or nearest_gene (aliases: nearest, ng)"
  )
}

normalize_vmr_set <- function(vmr_set = "shared") {
  x <- tolower(trimws(as.character(vmr_set)))
  x <- gsub("-", "_", x)
  dplyr::case_when(
    x %in% c("", "shared", "matched", "matched_shared", "current") ~ "shared",
    x %in% c("aa", "aa_only", "aa_all", "aa_specific", "population_aa") ~ "AA_only",
    x %in% c("ea", "ea_only", "ea_all", "ea_specific", "population_ea") ~ "EA_only",
    TRUE ~ NA_character_
  )
}

is_shared_plot_population <- function(population) {
  toupper(trimws(as.character(population))) %in% c("SHARED", "BOTH")
}

should_skip_shared_duplicate_population <- function(
    population,
    vmr_set = "shared",
    canonical = SHARED_VMR_CANONICAL_POPULATION) {
  normalize_vmr_set(vmr_set) == "shared" &&
    !is_shared_plot_population(population) &&
    toupper(trimws(as.character(population))) != toupper(canonical)
}

regctx_assoc_source_population <- function(
    population,
    vmr_set = "shared",
    canonical = SHARED_VMR_CANONICAL_POPULATION) {
  if (should_skip_shared_duplicate_population(population, vmr_set, canonical)) {
    return(toupper(canonical))
  }
  toupper(trimws(as.character(population)))
}

shared_vmr_plot_populations <- function(population, vmr_set = "shared") {
  vmr_set <- normalize_vmr_set(vmr_set)
  if (vmr_set != "shared") {
    return(toupper(trimws(as.character(population))))
  }
  if (is_shared_plot_population(population)) {
    return(c("AA", "EA"))
  }
  toupper(trimws(as.character(population)))
}

validate_vmr_set <- function(cohort, population, vmr_set = "shared") {
  vmr_set <- normalize_vmr_set(vmr_set)
  pop <- toupper(trimws(as.character(population)))
  if (is.na(vmr_set)) {
    stop("vmr_set must be shared, AA_only, or EA_only")
  }
  if (cohort == "BA_only" && vmr_set != "shared") {
    stop("BA_only supports vmr_set shared only")
  }
  if (vmr_set == "AA_only" && pop != "AA") {
    stop("vmr_set AA_only must be run with population AA")
  }
  if (vmr_set == "EA_only" && pop != "EA") {
    stop("vmr_set EA_only must be run with population EA")
  }
  if (vmr_set == "shared" && !pop %in% c("AA", "EA", "SHARED", "BOTH")) {
    stop("vmr_set shared requires population AA, EA, or SHARED (combined plots)")
  }
  vmr_set
}

regctx_output_dir <- function(cohort, tissue, population, modality = NULL,
                              run_tag = NULL, vmr_set = "shared") {
  vmr_set <- validate_vmr_set(cohort, population, vmr_set)
  parts <- c(
    here("heritability", "elastic_net_model", cohort,
         "tissue_comparison", "regulatory_context", "_m"),
    if (vmr_set == "shared") character(0) else vmr_set,
    tissue,
    population,
    modality,
    run_tag
  )
  do.call(file.path, as.list(parts[!is.na(parts) & nzchar(parts)]))
}

regctx_fig_dir <- function(cohort, population, vmr_set = "shared") {
  vmr_set <- validate_vmr_set(cohort, population, vmr_set)
  pop_dir <- if (vmr_set == "shared" && is_shared_plot_population(population)) {
    "shared"
  } else {
    toupper(trimws(as.character(population)))
  }
  parts <- c(
    here("heritability", "elastic_net_model", cohort,
         "tissue_comparison", "regulatory_context", "_m", "figures"),
    if (vmr_set == "shared") character(0) else vmr_set,
    pop_dir
  )
  do.call(file.path, as.list(parts[!is.na(parts) & nzchar(parts)]))
}

ensure_vmr_set_column <- function(df, vmr_set = "shared") {
  if (!nrow(df)) return(df)
  if (!"vmr_set" %in% names(df)) df$vmr_set <- "shared"
  df$vmr_set[is.na(df$vmr_set) | df$vmr_set == ""] <- "shared"
  df$vmr_set <- vapply(df$vmr_set, normalize_vmr_set, character(1))
  df
}

classify_h2 <- function(h2_unscaled, r_squared_cv) {
  factor(
    dplyr::case_when(
      r_squared_cv <= 0.3 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.3 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.3 ~ "Heritable",
      TRUE ~ NA_character_
    ),
    levels = H2_GROUP_LEVELS
  )
}

## Column names for model formulas: backtick-quote so names like D1-MSN or D2-SPN
## are single predictors (not subtraction / stats::D()).
quote_formula_sym <- function(nm) {
  sprintf("`%s`", gsub("`", "\\`", nm, fixed = TRUE))
}

## Row names (feature IDs) that have enough finite values, are not all ~zero,
## and are not constant — association tests cannot identify signal otherwise.
feature_ids_with_variation <- function(mat, sample_ids,
                                       zero_tol = 1e-12) {
  m <- mat[, sample_ids, drop = FALSE]
  if (nrow(m) == 0 || ncol(m) == 0) {
    return(character(0))
  }
  ok <- apply(m, 1L, function(x) {
    xf <- x[is.finite(x)]
    length(xf) >= 2L &&
      max(abs(xf)) >= zero_tol &&
      stats::var(xf) > 0
  })
  rownames(m)[ok]
}

normalize_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

safe_fwrite <- function(x, fn, ...) {
  dir.create(dirname(fn), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, fn, ...)
}

safe_read <- function(fn, ...) {
  if (!file.exists(fn)) stop("Missing required file: ", fn)
  data.table::fread(fn, ...)
}

pick_col <- function(df, candidates, required = TRUE, label = NULL) {
  nm <- intersect(candidates, colnames(df))
  if (length(nm) > 0) return(nm[[1]])
  if (required) {
    stop("Missing required column",
         if (!is.null(label)) paste0(" for ", label) else "",
         ". Tried: ", paste(candidates, collapse = ", "))
  }
  NA_character_
}

load_pheno <- function(cohort, tissue) {
  pheno_file <- if (cohort == "BA_only") {
    here("inputs", "phenotypes", "_m", "phenotypes-AA.tsv")
  } else {
    here("inputs", "phenotypes", "_m", "phenotypes-all.tsv")
  }

  pheno <- safe_read(pheno_file, na.strings = c("", "NA", "NaN"))
  pheno |>
    filter(region == tissue, primarydx == "Control", agedeath > 17) |>
    mutate(
      brnum = as.character(brnum),
      age = agedeath,
      dx = primarydx,
      education = case_when(
        education %in% c("7th", "8th", "Less than 7th",
                         "9th", "10th", "11th", "12th") ~ "less_than_hs",
        education %in% c("H.S. diploma", "GED") ~ "hs",
        education %in% c("1 yr college", "3 yrs college",
                         "Associate's or 2 yrs college", "Bachelor's",
                         "Master's", "JD", "PhD") ~ "more_than_hs",
        TRUE ~ as.character(education)
      ),
      marital_status = case_when(
        marital_status %in% c("Single") ~ "single",
        marital_status %in% c("Married") ~ "married",
        marital_status %in% c("Divorced", "Separated", "Widowed") ~
          "previously_married",
        TRUE ~ as.character(marital_status)
      ),
      any_trauma_hx = case_when(
        hx_sexual_abuse | hx_physical_abuse | hx_other_trauma |
          hx_military_service ~ 1L,
        is.na(hx_sexual_abuse) & is.na(hx_physical_abuse) &
          is.na(hx_other_trauma) & is.na(hx_military_service) ~ NA_integer_,
        TRUE ~ 0L
      )
    )
}

load_cell_props <- function(tissue) {
  fn <- here("inputs", "cell_proportions", "_m",
             paste0("music-proportions-", tissue, ".tsv"))
  prop <- safe_read(fn)
  prop |>
    mutate(sample_id = as.character(sample_id)) |>
    pivot_wider(names_from = "cell_type", values_from = "proportion") |>
    as.data.frame() |>
    tibble::column_to_rownames("sample_id")
}

load_rse <- function(tissue, modality) {
  fn <- here("inputs", "counts", rna_file_map[[modality]][[tissue]])
  if (!file.exists(fn)) stop("Missing RSE file: ", fn)
  env <- new.env(parent = emptyenv())
  load(fn, envir = env)
  objs <- mget(ls(env), envir = env)
  is_se <- vapply(objs, function(x) inherits(x, "SummarizedExperiment"), logical(1))
  if (!any(is_se)) stop("No SummarizedExperiment object found in ", fn)
  objs[[which(is_se)[1]]]
}

sample_ids_from_rse <- function(rse) {
  cd <- as.data.frame(colData(rse))
  if ("BrNum" %in% colnames(cd)) return(as.character(cd$BrNum))
  if ("brnum" %in% colnames(cd)) return(as.character(cd$brnum))
  as.character(colnames(rse))
}

get_assay_matrix <- function(rse, modality) {
  an <- assayNames(rse)
  assay_name <- if (modality == "expression") {
    dplyr::case_when(
      "logcounts" %in% an ~ "logcounts",
      "vst" %in% an ~ "vst",
      "tpm" %in% an ~ "tpm",
      "counts" %in% an ~ "counts",
      TRUE ~ an[[1]]
    )
  } else {
    dplyr::case_when(
      "psi" %in% an ~ "psi",
      "PSI" %in% an ~ "PSI",
      TRUE ~ an[[1]]
    )
  }
  mat <- as.matrix(assay(rse, assay_name))
  storage.mode(mat) <- "numeric"

  if (modality == "expression" && assay_name == "counts") {
    lib <- colSums(mat, na.rm = TRUE)
    lib[lib <= 0] <- NA_real_
    mat <- log2(t(t(mat) / lib * 1e6) + 1)
  }
  attr(mat, "assay_name") <- assay_name
  mat
}

feature_annotation <- function(rse, modality) {
  rd <- as.data.frame(rowData(rse))
  if (modality == "expression") {
    gene_id_col <- pick_col(rd, c("gene_id", "Geneid", "geneID"), label = "gene_id")
    gene_name_col <- pick_col(rd, c("gene_name", "symbol", "gene_symbol"),
                              required = FALSE, label = "gene_name")
    tibble(
      feature_id = rownames(rse),
      gene_id = as.character(rd[[gene_id_col]]),
      gene_id_base = strip_ensembl_version(rd[[gene_id_col]]),
      gene_name = if (!is.na(gene_name_col)) as.character(rd[[gene_name_col]]) else NA_character_
    )
  } else {
    tibble(
      feature_id = rownames(rse),
      psi_uid = rownames(rse)
    )
  }
}

metadata_for_residualization <- function(rse, pheno, cell_props) {
  cd <- as.data.frame(colData(rse))
  cd$sample_id <- sample_ids_from_rse(rse)
  cd$sample_id <- as.character(cd$sample_id)

  meta <- cd |>
    inner_join(pheno, by = c("sample_id" = "brnum"))
  cp <- cell_props |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id")
  meta <- meta |>
    left_join(cp, by = "sample_id")

  age_col <- pick_col(meta, c("Age", "age", "agedeath"), label = "Age")
  sex_col <- pick_col(meta, c("Sex", "sex"), label = "Sex")
  rin_col <- pick_col(meta, c("RIN", "rin", "RINscore", "RIN_sum"), label = "RIN")
  mod_col <- pick_col(meta, c("MoD", "mod", "manner_of_death"), label = "MoD")
  mito_col <- pick_col(meta, c("mito_mapping_rate", "mito_rate",
                               "mitochondrial_mapping_rate"), label = "mito_mapping_rate")
  assigned_col <- pick_col(meta, c("percent_assigned", "pct_assigned",
                                   "assigned_percent"), label = "percent_assigned")

  cell_cols <- setdiff(colnames(cp), "sample_id")
  if (length(cell_cols) == 0) stop("No cell-type proportion columns found.")

  out <- meta |>
    transmute(
      sample_id = sample_id,
      Age = as.numeric(.data[[age_col]]),
      Sex = as.factor(.data[[sex_col]]),
      RIN = as.numeric(.data[[rin_col]]),
      MoD = as.factor(.data[[mod_col]]),
      mito_mapping_rate = as.numeric(.data[[mito_col]]),
      percent_assigned = as.numeric(.data[[assigned_col]])
    )

  cells <- meta[, c("sample_id", cell_cols), drop = FALSE]
  for (cc in cell_cols) {
    cells[[cc]] <- asin(sqrt(pmax(0, pmin(1, as.numeric(cells[[cc]])))))
  }

  left_join(out, cells, by = "sample_id")
}

residualize_matrix <- function(mat, meta) {
  keep_samples <- intersect(colnames(mat), meta$sample_id)
  if (length(keep_samples) < 20) {
    stop("Too few overlapping samples for residualization: ", length(keep_samples))
  }
  mat <- mat[, keep_samples, drop = FALSE]
  meta <- meta[match(keep_samples, meta$sample_id), , drop = FALSE]

  covars <- setdiff(colnames(meta), "sample_id")
  complete <- complete.cases(meta[, covars, drop = FALSE])
  meta <- meta[complete, , drop = FALSE]
  mat <- mat[, complete, drop = FALSE]

  form <- as.formula(paste(
    "~",
    paste(vapply(covars, quote_formula_sym, character(1)), collapse = " + ")
  ))
  design <- model.matrix(form, data = meta)
  qr_mod <- qr(design)
  if (qr_mod$rank < ncol(design)) {
    stop("Residualization design is rank deficient: rank ", qr_mod$rank,
         " < ", ncol(design))
  }

  res <- t(apply(mat, 1, function(y) {
    ok <- is.finite(y)
    if (sum(ok) <= ncol(design) + 2) return(rep(NA_real_, length(y)))
    fit <- lm.fit(design[ok, , drop = FALSE], y[ok])
    out <- rep(NA_real_, length(y))
    out[ok] <- fit$residuals
    out
  }))
  colnames(res) <- meta$sample_id
  rownames(res) <- rownames(mat)
  list(residuals = res, meta = meta, design = design)
}

load_methylation <- function(cohort, tissue, pheno) {
  base <- if (cohort == "BA_only") {
    here("vmr-analysis", tissue, "_m", "vmr")
  } else {
    here("vmr-analysis", "all_individuals", tissue, "_m", "vmr")
  }
  files <- list.files(base, pattern = "_meth\\.phen$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop("No methylation phenotype files found under ", base)

  meth <- rbindlist(lapply(files, function(fn) {
    chr <- basename(dirname(fn))
    pos <- strsplit(sub("_meth\\.phen$", "", basename(fn)), "_")[[1]]
    df <- fread(fn, select = c(1, 3), header = FALSE)
    colnames(df) <- c("brnum", "meth")
    df |>
      mutate(
        brnum = as.character(brnum),
        seqnames = normalize_chr(sub("^chr_", "", chr)),
        start = as.integer(pos[[1]]),
        end = as.integer(pos[[2]]),
        vmr_id = paste(seqnames, start, end, sep = "_")
      )
  }), use.names = TRUE)

  meth |>
    inner_join(pheno |> dplyr::select(brnum), by = "brnum")
}

load_enet <- function(cohort, tissue, population = "AA", vmr_set = "shared") {
  vmr_set <- validate_vmr_set(cohort, population, vmr_set)
  if (cohort == "BA_only") {
    fn <- here("heritability", "elastic_net_model", "BA_only", tissue, "_m",
               paste0(tissue, "_summary_elastic-net.tsv"))
    safe_read(fn) |>
      mutate(
        seqnames = normalize_chr(chrom),
        h2_category = classify_h2(h2_unscaled, r_squared_cv),
        population = "AA",
        vmr_set = vmr_set
      ) |>
      dplyr::select(seqnames, start, end, h2_category, h2_unscaled,
                    r_squared_cv, num_snps, population, vmr_set)
  } else {
    if (vmr_set == "shared") {
      fn <- here("heritability", "elastic_net_model", "all_individuals",
                 tissue, "_m",
                 paste0(tissue, "_summary_elastic-net_matched_r2_0.3.tsv"))
      h2_col <- paste0("h2_unscaled_", population)
      r2_col <- paste0("r_squared_cv_", population)
      snp_col <- paste0("num_snps_", population)
      safe_read(fn) |>
        dplyr::rename(seqnames = chrom) |>
        mutate(
          h2_unscaled = .data[[h2_col]],
          r_squared_cv = .data[[r2_col]],
          num_snps = .data[[snp_col]],
          ## Matched-set h2_category is concordant across AA/EA (see tissue_comparison/h2_distribution).
          h2_category = factor(h2_category, levels = H2_GROUP_LEVELS),
          population = population,
          vmr_set = vmr_set
        ) |>
        dplyr::select(seqnames, start, end, h2_category, h2_unscaled,
                      r_squared_cv, num_snps, population, vmr_set)
    } else {
      fn <- here("heritability", "elastic_net_model", "all_individuals",
                 tissue, "_m",
                 paste0(tissue, "_summary_elastic-net_", population, ".tsv"))
      safe_read(fn) |>
        dplyr::rename(seqnames = chrom) |>
        mutate(
          seqnames = normalize_chr(seqnames),
          h2_category = classify_h2(h2_unscaled, r_squared_cv),
          population = population,
          vmr_set = vmr_set
        ) |>
        dplyr::select(seqnames, start, end, h2_category, h2_unscaled,
                      r_squared_cv, num_snps, population, vmr_set)
    }
  }
}

load_abc_links <- function(cohort, population = "AA") {
  fn <- if (cohort == "BA_only") {
    here("heritability", "elastic_net_model", "BA_only", "tissue_comparison",
         "annotation", "open_chromatin", "_m", "abc_vmr_gene_links.tsv")
  } else {
    here("heritability", "elastic_net_model", "all_individuals", "tissue_comparison",
         "annotation", "open_chromatin", "_m",
         paste0("abc_vmr_gene_links_", population, ".tsv"))
  }
  safe_read(fn) |>
    filter(tissue != "Pooled") |>
    mutate(
      tissue = tolower(tissue),
      seqnames = normalize_chr(seqnames),
      TargetGene_base = strip_ensembl_version(TargetGene),
      vmr_id = paste(seqnames, start, end, sep = "_")
    )
}

load_nearest_gene_links <- function(cohort, tissue, population = "AA",
                                    window = 250000L,
                                    vmr_set = "shared") {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install org.Hs.eg.db for nearest-gene links (Bioconductor).")
  }
  suppressPackageStartupMessages(library("org.Hs.eg.db", character.only = TRUE))
  annot_dir <- here(
    "heritability", "elastic_net_model", cohort,
    "tissue_comparison", "annotation", "_m"
  )
  fn <- file.path(
    annot_dir,
    paste0(tissue, "_vmr_genes_within_250kb_hg38.tsv")
  )
  df <- safe_read(fn)

  sn_col <- intersect(c("seqnames", "chr", "chrom"), colnames(df))
  if (length(sn_col) == 0) {
    stop("No seqnames/chr column in nearest-gene annotation: ", fn)
  }
  df$seqnames <- normalize_chr(df[[sn_col[[1]]]])
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)

  eg_col <- intersect(
    c("nearest_gene_id_within_250kb", "nearest_gene_id"),
    colnames(df)
  )
  if (length(eg_col) == 0) {
    stop("Missing nearest Entrez column (expected nearest_gene_id_within_250kb): ", fn)
  }
  eg_col <- eg_col[[1]]

  sym_col <- intersect(
    c("nearest_gene_symbol_within_250kb", "nearest_gene_symbol"),
    colnames(df)
  )
  sym_col <- if (length(sym_col)) sym_col[[1]] else NA_character_

  dist_col <- intersect(
    c(
      "distance_to_nearest_gene_within_250kb",
      "distance_to_nearest_gene"
    ),
    colnames(df)
  )
  dist_col <- if (length(dist_col)) dist_col[[1]] else NA_character_

  df$vmr_id <- paste(df$seqnames, df$start, df$end, sep = "_")
  df$nearest_entrez <- as.character(df[[eg_col]])
  df <- df[!is.na(df$nearest_entrez) & df$nearest_entrez != "", , drop = FALSE]

  if (nrow(df) == 0) return(tibble())

  eg_unique <- unique(df$nearest_entrez)
  map_ens <- suppressWarnings(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = eg_unique,
      column = "ENSEMBL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
  )
  entrez_to_base <- strip_ensembl_version(as.character(map_ens))
  names(entrez_to_base) <- names(map_ens)

  df$gene_id_base <- entrez_to_base[df$nearest_entrez]
  df <- df[!is.na(df$gene_id_base) & df$gene_id_base != "", , drop = FALSE]

  if (!is.na(dist_col)) {
    d <- as.numeric(df[[dist_col]])
    df <- df[is.na(d) | d <= as.numeric(window), , drop = FALSE]
  }

  enet <- load_enet(cohort, tissue, population, vmr_set) |>
    mutate(vmr_id = paste(seqnames, start, end, sep = "_"))
  df <- df |>
    dplyr::select(-any_of(c("h2_category", "h2_unscaled"))) |>
    tibble::as_tibble() |>
    inner_join(
      enet |> dplyr::select(vmr_id, h2_category, h2_unscaled),
      by = "vmr_id"
    )

  out <- tibble::tibble(
    seqnames = df$seqnames,
    start = df$start,
    end = df$end,
    vmr_id = df$vmr_id,
    h2_category = df$h2_category,
    h2_unscaled = as.numeric(df$h2_unscaled),
    gene_id_base = df$gene_id_base,
    nearest_gene_symbol_within_250kb = if (!is.na(sym_col)) {
      as.character(df[[sym_col]])
    } else {
      NA_character_
    },
    distance_to_nearest_gene_within_250kb = if (!is.na(dist_col)) {
      as.numeric(df[[dist_col]])
    } else {
      NA_real_
    },
    tissue = tissue,
    population = population,
    vmr_set = vmr_set
  )

  distinct(out)
}

load_psi_links <- function(cohort, tissue, population = "AA", window = 250000,
                           vmr_set = "shared") {
  enet <- load_enet(cohort, tissue, population, vmr_set) |>
    mutate(vmr_id = paste(seqnames, start, end, sep = "_"))
  psi_annot <- safe_read(here("inputs", "counts", "psi-annotation.tsv"))
  if (!"psi_uid" %in% colnames(psi_annot)) {
    psi_annot$psi_uid <- rownames(psi_annot)
  }
  psi_annot <- psi_annot |>
    mutate(
      chrom = normalize_chr(chrom),
      psi_feature = psi_uid
    )

  vmr_gr <- GRanges(enet$seqnames, IRanges(enet$start, enet$end))
  psi_gr <- GRanges(psi_annot$chrom,
                    IRanges(pmax(1L, psi_annot$start - window),
                            psi_annot$end + window))
  hits <- findOverlaps(vmr_gr, psi_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(tibble())

  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  psi_event_gr <- GRanges(
    psi_annot$chrom[sh],
    IRanges(psi_annot$start[sh], psi_annot$end[sh])
  )
  tibble(
    seqnames = enet$seqnames[qh],
    start = enet$start[qh],
    end = enet$end[qh],
    vmr_id = enet$vmr_id[qh],
    h2_category = enet$h2_category[qh],
    h2_unscaled = enet$h2_unscaled[qh],
    r_squared_cv = enet$r_squared_cv[qh],
    num_snps = enet$num_snps[qh],
    tissue = tissue,
    population = population,
    vmr_set = vmr_set,
    psi_uid = psi_annot$psi_feature[sh],
    gene_id = psi_annot$gene_id[sh],
    gene_name = psi_annot$gene_name[sh],
    event_type = psi_annot$event_type[sh],
    distance = pmax(0, distance(vmr_gr[qh], psi_event_gr))
  )
}

load_architecture_covariates <- function(cohort, tissue, population = "AA",
                                         vmr_set = "shared") {
  annot_dir <- here("heritability", "elastic_net_model", cohort,
                    "tissue_comparison", "annotation", "_m")
  annot <- safe_read(file.path(annot_dir,
                               paste0(tissue, "_vmr_annotations_hg38_wide.tsv"))) |>
    mutate(
      seqnames = normalize_chr(seqnames),
      vmr_length = end - start,
      vmr_id = paste(seqnames, start, end, sep = "_")
    ) |>
    dplyr::select(-any_of(c(
      "h2_category", "h2_unscaled", "r_squared_cv", "num_snps", "population"
    )))
  enet <- load_enet(cohort, tissue, population, vmr_set) |>
    mutate(vmr_id = paste(seqnames, start, end, sep = "_"))

  out <- annot |>
    left_join(enet, by = c("seqnames", "start", "end", "vmr_id")) |>
    mutate(
      chromatin_state = case_when(
        hg38_genes_promoters == 1 ~ "Promoter",
        hg38_enhancers_fantom == 1 ~ "Enhancer",
        hg38_genes_intergenic == 1 ~ "Intergenic",
        TRUE ~ "Genic/other"
      )
    )
  cpg_cols <- intersect(c("hg38_cpg_islands", "hg38_cpg_shores",
                          "hg38_cpg_shelves"), colnames(out))
  out$cpg_density_proxy <- if (length(cpg_cols) > 0) {
    rowSums(as.data.frame(out)[, cpg_cols, drop = FALSE], na.rm = TRUE)
  } else {
    NA_real_
  }
  out
}
