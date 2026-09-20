#!/usr/bin/env Rscript
#### 08_region_donor_generalization -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA [--allow-unlocked]
##
## A run of this module is NOT per-region: the deliverable IS the comparison
## across regions, so one run holds all three and `region` is the literal
## "crossregion". Splitting it per region would make the primary analysis a
## thing assembled outside any run directory, which is the pattern AGENTS.md 9
## exists to prevent.
##
## Every upstream this module reads is gated here, at creation time, and its run
## ID is written into the manifest. Nothing downstream re-checks acceptance for
## a sealed run (that is require_accepted_upstream's contract), so this is the
## only place the dependency graph is proven.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "08_region_donor_generalization"
MODULE_TAG <- "rdg"

opts <- parse_v2_args(require = c("cohort"))
allow_unlocked <- isTRUE(opts$allow_unlocked)
cohort <- opts$cohort

cfg <- load_config("region_donor_generalization")
thresholds <- load_config("analysis_thresholds")
assert_locked(list(region_donor_generalization = cfg),
              allow_unlocked = allow_unlocked)

## The region axis runs in a discovery ARM. A cell would be wrong here: its
## loci are a different set from the arm's, and 04-07 are accepted on the arm.
if (!identical(parse_cell(cohort)$cell_kind, "arm")) {
    stop("The region axis runs in a discovery arm, not the estimation cell '",
         cohort, "'. The donor-group axis reads its cells from ",
         "_m/combined/ and does not need a cell token here.")
}
regions <- as.character(config_get(cfg, "cross_region_replication.regions"))

## ------------------------------------------------- region axis: gate 01-07
##
## AGENTS.md 6 names 04 and 05 as this module's blocking upstreams. 02, 03 and
## 07 are gated too, because tier 1 reads all three: the replication claim is
## about the genetic-control architecture, and 02's score, 03's prediction and
## 07's coupling are three of its four legs.
region_upstreams <- c("01_vmr_catalog", "02_local_genetic_variance",
                      "03_local_snp_prediction",
                      "04_repeat_repressive_architecture",
                      "05_cpg_meqtl_burden",
                      "07_transcription_splicing_coupling")

accepted <- list()
for (mod in region_upstreams) {
    for (re in regions) {
        row <- require_accepted_upstream(mod, cohort, re,
                                         allow_unaccepted = allow_unlocked)
        accepted[[paste(mod, re, sep = "|")]] <- row
    }
}

## One locus set across the whole region axis, or "replication across regions"
## is comparing different loci that happen to share an ID scheme. Each region
## has its OWN vmr_set_id by construction (discovery is per region), so the
## assertion is within region, across modules.
for (re in regions) {
    ids <- unique(na.omit(vapply(
        region_upstreams,
        function(m) as.character(accepted[[paste(m, re, sep = "|")]]$vmr_set_id),
        character(1))))
    if (length(ids) > 1L) {
        stop("Upstream modules cite different vmr_set_ids for ", cohort, " x ",
             re, ": ", paste(ids, collapse = ", "),
             ". Cross-region replication cannot be assembled from ",
             "inconsistent locus sets.")
    }
}

## ------------------------------------------- donor-group axis: gate the cells
##
## The cells are gated, and so is the policy. A widened policy must stop the
## run here rather than at the figure stage.
policy <- donor_group_inference_policy()
dg_cells <- as.character(config_get(cfg, "donor_group.cells"))
dg_regions <- as.character(config_get(cfg, "donor_group.regions"))

cell_runs <- list()
for (cell in dg_cells) {
    for (re in dg_regions) {
        for (mod in c("01b_estimation_cells", "02_local_genetic_variance",
                      "03_local_snp_prediction")) {
            row <- require_accepted_upstream(mod, cell, re,
                                             allow_unaccepted = allow_unlocked)
            cell_runs[[paste(mod, cell, re, sep = "|")]] <- row
        }
    }
}

## The design in one assertion: the two cells were discovered ONCE, so within a
## region they must carry the same vmr_set_id, and it must differ from the arm's.
for (re in dg_regions) {
    ids <- unique(unlist(lapply(dg_cells, function(cell) {
        as.character(cell_runs[[paste("02_local_genetic_variance", cell, re,
                                      sep = "|")]]$vmr_set_id)
    })))
    if (length(ids) != 1L) {
        stop("The donor-group cells carry different vmr_set_ids for ", re,
             ": ", paste(ids, collapse = ", "), " (AGENTS.md 7.7)")
    }
    arm_id <- as.character(
        accepted[[paste("02_local_genetic_variance", re, sep = "|")]]$vmr_set_id)
    if (identical(ids, arm_id)) {
        stop("The donor-group cells share the '", cohort, "' arm's vmr_set_id ",
             "for ", re, ". They must sit on the pooled catalog, or the axis ",
             "is the prohibited AA-vs-all_individuals contrast.")
    }
}

## The recombination outputs must already exist: 08 reads them, it does not
## recompute a score (AGENTS.md 7.6 forbids re-ranking across cells).
dg_inputs <- config_get(cfg, "donor_group.inputs")
for (nm in names(dg_inputs)) {
    for (re in dg_regions) {
        f <- file.path(V2_ROOT, gsub("{region}", re, dg_inputs[[nm]],
                                     fixed = TRUE))
        if (!file.exists(f)) {
            stop("Donor-group input '", nm, "' is missing for ", re, ": ", f,
                 "\n  Run 02/_h/14_combine_donor_group_cells.R and ",
                 "03/_h/07_stack_donor_group_predictions.R first.")
        }
    }
}

## ------------------------------------------ tier 3: gate the subsample cells
##
## Tier 3 is the one axis that needs runs this module's own config commissioned.
## Gate them exactly like any other upstream: a downsampling conclusion drawn
## from an unaccepted run is not citable (AGENTS.md 6).
ds <- config_get(cfg, "caudate_downsampling")
ds_runs <- list()
if (isTRUE(ds$enabled)) {
    ds_region <- as.character(ds$region)
    ## cell_token_template is the FULL token, arm included, so substitute all
    ## three keys and prepend nothing -- prepending source_cell as well would
    ## build 'AA.AA.n118r1'.
    ds_cells <- vapply(seq_len(as.integer(ds$n_replicates)), function(i) {
        tok <- ds$cell_token_template
        tok <- gsub("{source_cell}", ds$source_cell, tok, fixed = TRUE)
        tok <- gsub("{target_n}", ds$target_n, tok, fixed = TRUE)
        gsub("{replicate}", i, tok, fixed = TRUE)
    }, character(1))
    for (cell in ds_cells) {
        ## Fail early and legibly if the cell is not even declared, rather than
        ## with an opaque parse error deeper in.
        parsed <- parse_cell(cell)
        if (!identical(parsed$cell_kind, "donor_subsample")) {
            stop("Tier-3 cell '", cell, "' is not a donor_subsample cell")
        }
        drawn_n <- as.integer(sub("^n(\\d+)r\\d+$", "\\1",
                                  parsed$estimation_group))
        if (!identical(drawn_n, as.integer(ds$target_n))) {
            stop("Tier-3 cell '", cell, "' draws n=", drawn_n, ", not ",
                 ds$target_n)
        }
        for (mod in c("01b_estimation_cells", "02_local_genetic_variance",
                      "03_local_snp_prediction")) {
            row <- require_accepted_upstream(mod, cell, ds_region,
                                             allow_unaccepted = allow_unlocked)
            ## Key on mod|cell|REGION, exactly like `cell_runs` above. Keying on
            ## mod|cell alone wrote `upstream_03_local_snp_prediction_AA_n118r1`
            ## while Stage 04 looks up `..._AA_n118r1_caudate`, so every
            ## replicate missed and fell through Stage 04's arm fallback to the
            ## full caudate run -- comparing the arm against itself three times.
            ds_runs[[paste(mod, cell, ds_region, sep = "|")]] <- row
        }
    }
    ## The subsets must sit on the SAME locus set as the full caudate run, or
    ## the comparison confounds donor count with locus turnover -- which would
    ## defeat the entire point of the tier.
    full_id <- as.character(
        accepted[[paste("02_local_genetic_variance", ds_region,
                        sep = "|")]]$vmr_set_id)
    for (cell in ds_cells) {
        sub_id <- as.character(
            ds_runs[[paste("02_local_genetic_variance", cell, ds_region,
                           sep = "|")]]$vmr_set_id)
        if (!identical(sub_id, full_id)) {
            stop("Tier-3 cell '", cell, "' sits on vmr_set_id ", sub_id,
                 " but the full ", ds_region, " run sits on ", full_id,
                 ". A donor-count sensitivity must hold the locus set fixed.")
        }
    }
}

## ------------------------------------------------------------ open the run
flatten <- function(lst, prefix) {
    if (!length(lst)) return(list())
    stats::setNames(
        lapply(lst, function(x) as.character(x$run_id %||% NA_character_)),
        paste0(prefix, gsub("[|.]", "_", names(lst))))
}

run <- new_run(
    module = MODULE_TAG,
    cohort = cohort,
    region = "crossregion",
    module_root = file.path(V2_ROOT, MODULE),
    ## Deliberately NA: a run spanning three regions has three locus sets, and
    ## writing one of them here would misattribute the other two. They are
    ## recorded per region in the manifest fields below.
    vmr_set_id = NA_character_,
    extra = c(
        list(
            smoke_run          = allow_unlocked,
            regions            = paste(regions, collapse = ","),
            donor_group_cells  = paste(dg_cells, collapse = ","),
            donor_group_inference = policy$mode,
            ancestry_effect_claim_allowed = policy$ancestry_effect_claim_allowed,
            cross_group_raw_score_comparison =
                policy$cross_group_raw_score_comparison,
            tier3_enabled      = isTRUE(ds$enabled),
            tier3_target_n     = as.integer(ds$target_n),
            tier3_replicates   = as.integer(ds$n_replicates),
            config_region_donor_generalization_sha256 =
                attr(cfg, "config_sha256"),
            config_analysis_thresholds_sha256 = policy$config_sha256
        ),
        stats::setNames(
            lapply(regions, function(re) as.character(
                accepted[[paste("02_local_genetic_variance", re,
                                sep = "|")]]$vmr_set_id)),
            paste0("vmr_set_id_", regions)),
        flatten(accepted, "upstream_"),
        flatten(cell_runs, "upstream_"),
        flatten(ds_runs, "upstream_")
    ))

dir.create(file.path(run$dir, "results"), recursive = TRUE,
           showWarnings = FALSE)

## The tier table travels with the run, so a reader of the run directory alone
## knows what each output is licensed to say.
tier_tbl <- rbindlist(lapply(names(cfg$tiers), function(nm) {
    t <- cfg$tiers[[nm]]
    data.table(tier = nm, rank = t$rank, label = t$label,
               licenses = trimws(t$licenses),
               regions = paste(t$regions %||% "all", collapse = ","),
               cannot_establish = paste(t$cannot_establish %||% "",
                                        collapse = "; "))
}))[order(rank)]
write_atomic(tier_tbl, file.path(run$dir, "results", "tiers.tsv"))

message("[run] ", run$run_id, ": ", length(accepted), " region-axis + ",
        length(cell_runs), " cell + ", length(ds_runs),
        " tier-3 accepted upstreams")

## The launcher reads this line.
cat(run$run_id, "\n", sep = "")
