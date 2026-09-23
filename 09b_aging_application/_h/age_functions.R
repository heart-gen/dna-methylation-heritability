#### 09b_aging_application -- module-local functions ####
##
## Sourced by 01, 02 and 05. Everything here is linear algebra on a donor x VMR
## methylation matrix, because the donor bootstraps refit every VMR a thousand
## times: one QR per design, applied to all VMRs at once,
## replaces thousands of lm() calls with a few matrix products. The smoke test
## checks that these agree with per-VMR lm() to 1e-10.

## Read one Module 01 phenotype per eligible VMR into a donor x VMR matrix.
##
## Covariates come from 00_shared/locus_io.R::load_locus_phenotype(), the reader
## Module 02 used. Module 01 writes one covariate file per chromosome; they are
## asserted identical to the first chromosome's, so a single donor table is
## correct for every VMR. Each phenotype is ALIGNED to the reference donor order
## by ID, never assumed to share it (AGENTS.md 7.1).
load_age_inputs <- function(tasks, vmr_run_dir, catalog_cohort, covar_prefix) {
    ref <- NULL
    cov_seen <- list()
    Y <- matrix(NA_real_, nrow = 0L, ncol = nrow(tasks))
    for (i in seq_len(nrow(tasks))) {
        task <- as.list(tasks[i])
        ph <- load_locus_phenotype(task = task, vmr_run_dir = vmr_run_dir,
                                   cohort = catalog_cohort,
                                   covar_prefix = covar_prefix)
        if (!is.null(ph$pcs)) {
            stop("Estimation-cell genotype PCs present for ", task$vmr_id,
                 "; this module is prespecified for the AA arm only.")
        }
        chrom <- as.character(task$chrom)
        if (is.null(ref)) {
            ref <- merge(ph$covar, ph$qcovar, by = c("FID", "IID"), all = FALSE)
            ref <- ref[order(ref$FID), , drop = FALSE]
            if (anyDuplicated(ref$FID)) stop("Duplicate donors in Module 01 covariates")
            if (nrow(ref) != nrow(ph$covar) || nrow(ref) != nrow(ph$qcovar)) {
                stop("covar and qcovar donor sets differ for ", chrom)
            }
            Y <- matrix(NA_real_, nrow = nrow(ref), ncol = nrow(tasks),
                        dimnames = list(ref$FID, tasks$vmr_id))
        }
        if (is.null(cov_seen[[chrom]])) {
            this <- merge(ph$covar, ph$qcovar, by = c("FID", "IID"), all = FALSE)
            this <- this[order(this$FID), , drop = FALSE]
            rownames(this) <- NULL
            chk <- ref; rownames(chk) <- NULL
            if (!isTRUE(all.equal(this, chk, check.attributes = FALSE))) {
                stop("Module 01 covariates differ between chromosomes (", chrom,
                     "); one donor table cannot serve every VMR.")
            }
            cov_seen[[chrom]] <- TRUE
        }
        if (anyDuplicated(ph$phenotype$FID)) {
            stop("Duplicate donors in phenotype for ", task$vmr_id)
        }
        idx <- match(ref$FID, ph$phenotype$FID)
        if (anyNA(idx) || nrow(ph$phenotype) != nrow(ref)) {
            stop("Phenotype donors for ", task$vmr_id, " differ from the ",
                 "covariate donors (AGENTS.md 10.1: missing donors fail loudly).")
        }
        Y[, i] <- as.numeric(ph$phenotype$phenotype[idx])
    }
    list(Y = Y, donors = data.table::as.data.table(ref))
}

## Donor-level design for one spec. Returns the rows (donor indices into the
## full table), the model matrix, and the permutation strata. Age is in decades.
build_spec_design <- function(spec, donors, cfg, cellpcs = list()) {
    ac <- cfg$age_model
    unit <- as.numeric(ac$age_unit_years)
    ctrl <- as.character(ac$diagnosis_control_label)
    keep <- rep(TRUE, nrow(donors))
    if (identical(spec$donors, "controls")) keep <- keep & donors$diagnosis == ctrl
    if (!is.null(spec$min_age)) keep <- keep & donors$age >= as.numeric(spec$min_age)

    d <- data.frame(age = donors$age / unit,
                    sex = factor(donors$sex),
                    diagnosis = factor(donors$diagnosis))
    extra <- as.character(unlist(spec$extra_terms))
    for (term in extra) {
        src <- sub("_cellPC[0-9]+$", "", term)
        pcs <- cellpcs[[src]]
        if (is.null(pcs)) stop("Spec needs ", src, " cell PCs, which were not built")
        col <- sub("^.*_", "", term)
        col <- sub("cellPC", "PC", col)
        if (!col %in% colnames(pcs)) stop("No ", col, " in ", src, " cell PCs")
        v <- pcs[match(donors$FID, rownames(pcs)), col]
        d[[term]] <- v
        keep <- keep & is.finite(v)
    }
    terms <- setdiff(c("age", "sex", "diagnosis"), unlist(spec$drop_terms))
    terms <- c(terms, extra)
    d <- d[keep, , drop = FALSE]
    for (f in c("sex", "diagnosis")) if (f %in% terms) d[[f]] <- droplevels(d[[f]])
    X <- stats::model.matrix(stats::reformulate(terms), data = d)
    if (qr(X)$rank < ncol(X)) {
        stop("Design for spec is rank-deficient: ", paste(colnames(X), collapse = ","))
    }
    list(rows = which(keep), X = X, age_col = which(colnames(X) == "age"),
         strata = as.character(donors$diagnosis[keep]),
         n = sum(keep),
         n_controls = sum(donors$diagnosis[keep] == ctrl),
         n_cases = sum(donors$diagnosis[keep] != ctrl))
}

## Age coefficient, SE, t and p for every VMR at once. Y is n x V.
## The generic matrix fitter lives in 00_shared/axis_inference.R, which Module 10
## also uses; this keeps 09b's call sites and its age-specific naming.
fit_age_matrix <- function(X, Y, age_col) fit_scalar_matrix(X, Y, age_col)

## resample_rows(), prepare_axis(), axis_estimate(), block_jackknife_se()
## and combined_inference() moved to 00_shared/axis_inference.R on
## 2026-09-19, arithmetic unchanged, when Module 10 needed the same
## machinery (AGENTS.md 5.3: refactor shared logic rather than duplicate).
## 09b/tests/ checks the shared versions against the definitions this file
## carried at commit abb0c6789.

## Per-VMR outcomes from one set of age fits.
##
## debiased_sq_beta = beta_hat^2 - SE^2 is UNBIASED for beta^2 whatever the SE.
## That is the property the design rests on: a VMR with strong local SNP control
## carries its genetic variance in the age model's residual, so its SE is
## larger, and any outcome whose expectation depends on SE (|beta_hat|, its
## rank, |t|, -log10 p) couples to the score mechanically. See config/aging.yml.
age_outcomes <- function(fit) {
    list(debiased_sq_beta = fit$beta^2 - fit$se^2,
         signed_beta_age = fit$beta)
}


## scMD integration gate: `scmd_gate_passes()` was defined here and MOVED to
## 00_shared/cell_composition.R on 2026-09-23, unchanged, when Module 04 needed
## the identical decision for its own composition arms. 00_shared/load.R sources
## it, so callers in this module are unaffected; the point of the move is that
## there is one definition of the gate rather than two that can drift
## (AGENTS.md 5.3).


## Which Module 04 feature columns are computed FROM the DNAm scMD donor
## proportions, and therefore inherit scMD's integration gate.
##
## READ FROM THE FEATURE TABLE, never inferred from a column name. The name
## `cell_composition_r2` is fixed by PI-locked configuration
## (config/aging.yml:axis.arms, config/gwas_negative_controls.yml) while the
## quantity behind it is not: before 2026-09-23 Module 04 built it from
## dnam-scmd-proportions-{region}.tsv in every region; since the MuSiC
## correction it builds it from the RNA MuSiC proportions in every region and
## records the modality in `cell_composition_r2_source`, with the scMD value
## kept beside it in `cell_composition_r2_scmd` (NA where the gate fails). A
## constant here would have been right against one table and wrong against the
## other, and would have put "scMD" in the audit trail of a run whose covariate
## was MuSiC -- worse than the original defect, because it would be a false
## reason for a defensible action.
##
## WHAT THE GATE IS FOR. AGENTS.md 7.4 is asymmetric: RNA MuSiC adjustment is
## required in every region, DNAm scMD adjustment only "when the integration
## gate passes", which is caudate alone (total-neuron rho 0.72 caudate, -0.03
## dlpfc, -0.10 hippocampus). So a MuSiC-derived covariate gates everywhere and
## an scMD-derived one gates in caudate only. The selection is on the recorded
## modality; it was never really about the name.
##
## WHY A GATED ARM RATHER THAN A METHYLATION PROPERTY. The column is a derived
## scalar, so it can be read as a property of the VMR's methylation -- AGENTS.md
## 7.4 does ask for "cell-composition-associated methylation properties" as an
## adjustment. But the adjective is load-bearing. Where the proportions behind it
## do not track composition, the R^2 measures methylation variance shared with
## three PCs of a quantity that is not composition. As a GATING member its output
## is a boolean that moves the module's verdict, and the only claim that boolean
## can license is "the gradient is not cell composition". Where that is
## unavailable the arm must be treated exactly as `cell_scmd` is: absent
## evidence, not contrary evidence.
##
## The gate is still applied at the point of CONSUMPTION, not construction. The
## column is well defined in every region whichever modality built it; what is
## region-conditional is whether it can support a composition claim.
##
## The DESCRIPTIVE annotation block (config/aging.yml:annotation_associations)
## keeps the column in every region on purpose, with the caveat written into the
## config at line 211: there it is the annotation being described, not an
## adjustment claiming to have removed composition, and it never reaches the
## region reading.

## Recognised values of `cell_composition_r2_source`. An unrecognised value is
## an error, not a default: this module would be deciding a gate over a modality
## it has never been told about.
SCMD_MODALITY <- c(rna_music = FALSE, dnam_scmd = TRUE)

## Used ONLY for a Module 04 table that predates `cell_composition_r2_source`
## (runs up to rra-AA-*-20260906). Such a table genuinely is scMD-derived, so the
## fallback is the correct answer for it -- but it is an assumption about data
## that cannot speak for itself, so it is announced rather than taken quietly.
SCMD_DERIVED_FEATURES_PRE_MODALITY <- c("cell_composition_r2")

#' @param feat Module 04's vmr-features.tsv, as a data.table.
#' @return the feature column names that are scMD-derived, possibly none.
scmd_derived_feature_columns <- function(feat, announce = message) {
    ## An explicitly named modality column means what it says whatever the
    ## table-level source is; both exist side by side after the correction.
    explicit <- intersect("cell_composition_r2_scmd", names(feat))
    if (!"cell_composition_r2_source" %in% names(feat)) {
        announce("[09b] Module 04 feature table carries no ",
                 "cell_composition_r2_source column: treating it as a ",
                 "pre-2026-09-23 table, in which cell_composition_r2 is ",
                 "DNAm scMD-derived. Provenance is ASSUMED, not read.")
        return(unique(c(SCMD_DERIVED_FEATURES_PRE_MODALITY, explicit)))
    }
    src <- unique(as.character(feat$cell_composition_r2_source))
    src <- src[!is.na(src)]
    if (!length(src)) {
        stop("cell_composition_r2_source is present but entirely NA; the ",
             "modality of cell_composition_r2 cannot be established")
    }
    if (length(src) > 1L) {
        stop("cell_composition_r2_source takes more than one value in one ",
             "feature table (", paste(src, collapse = ", "), "); one column ",
             "cannot have two provenances")
    }
    if (!src %in% names(SCMD_MODALITY)) {
        stop("Module 04 records cell_composition_r2_source '", src, "', which ",
             "this module does not recognise. Teach age_functions.R:",
             "SCMD_MODALITY whether that modality is scMD-derived; it is not ",
             "assumed either way.")
    }
    unique(c(if (isTRUE(SCMD_MODALITY[[src]])) "cell_composition_r2", explicit))
}

## Does an axis arm's covariate set or row subset depend on an scMD-derived
## feature? Arms are declared in config/aging.yml:axis.arms; `scmd_cols` comes
## from scmd_derived_feature_columns().
arm_is_scmd_derived <- function(arm, scmd_cols) {
    cols <- unique(c(as.character(unlist(arm$add_covariates)),
                     as.character(unlist(arm$subset$column))))
    any(cols %in% scmd_cols)
}

## An scMD-derived arm in the GATING role is not fitted where the scMD
## integration gate fails, mirroring the `requires_scmd_integration_gate` spec
## rule in 01_age_effects.R. A non-gating arm is still fitted: it cannot move
## the verdict, and it is flagged in axis-tests.tsv instead.
arm_skipped_for_scmd <- function(arm, scmd_ok, scmd_cols) {
    !isTRUE(scmd_ok) && identical(as.character(arm$role), "gating") &&
        arm_is_scmd_derived(arm, scmd_cols)
}


## Region reading under strict conjunction (config/aging.yml:region_reading).
##
## Extracted from 03_apply_gates.R, arithmetic unchanged, so tests can drive it
## from a sealed run's axis-tests.tsv without rerunning a stage -- and so the
## "conjunction over FITTED members only" rule has one implementation to check.
##
## `main` is axis-tests.tsv restricted to the primary outcome. `reason_for` maps
## a member name to why it was not fitted, for the audit trail.
region_reading_members <- function(main, rr, hyp,
                                   reason_for = function(m) NA_character_) {
    row_for <- function(member) {
        r <- main[(spec == member & arm == "base") |
                  (spec == "primary" & arm == member)]
        if (nrow(r) > 1L) stop("Ambiguous region-reading member: ", member)
        r
    }
    prim <- row_for("primary")
    if (!nrow(prim)) stop("The primary spec has no axis row")
    primary_supported <- prim$direction == hyp &&
        prim$p < as.numeric(rr$primary_alpha)

    member_rows <- list()
    for (m in as.character(unlist(rr$same_n_members))) {
        r <- row_for(m)
        if (!nrow(r)) {
            ## Not fitted in this region. Recorded with its reason, and it does
            ## not count against the conjunction: a member that cannot be
            ## fitted is absent evidence, not contrary evidence. Both reasons
            ## that reach here are the scMD integration gate -- for the
            ## `cell_scmd` spec (01_age_effects.R) and for any scMD-derived
            ## gating arm (02_axis_test.R).
            member_rows[[m]] <- data.table::data.table(
                member = m, rule = rr$same_n_rule, fitted = FALSE,
                survives = NA, estimate = NA_real_, p = NA_real_,
                reason = reason_for(m))
            next
        }
        member_rows[[m]] <- data.table::data.table(
            member = m, rule = rr$same_n_rule, fitted = TRUE,
            survives = r$direction == hyp && r$p < as.numeric(rr$primary_alpha),
            estimate = r$estimate, p = r$p, reason = NA_character_)
    }
    for (m in as.character(unlist(rr$reduced_n_members))) {
        r <- row_for(m)
        if (!nrow(r)) stop("Reduced-n member ", m, " was not fitted")
        frac <- r$estimate / prim$estimate
        member_rows[[m]] <- data.table::data.table(
            member = m, rule = rr$reduced_n_rule, fitted = TRUE,
            survives = r$direction == hyp && is.finite(frac) &&
                frac >= as.numeric(rr$reduced_n_min_fraction),
            estimate = r$estimate, p = r$p, reason = NA_character_,
            fraction_of_primary = frac)
    }
    members <- data.table::rbindlist(member_rows, fill = TRUE)

    ## The conjunction is over FITTED members only, and `all(logical(0))` is
    ## TRUE, so an empty same-n set would silently return a vacuous "survives
    ## every sensitivity". Refuse it rather than report it.
    same_n <- as.character(unlist(rr$same_n_members))
    if (!any(members$fitted & members$member %in% same_n)) {
        stop("No same-n gating member was fitted; the strict conjunction would ",
             "be vacuous. Check the scMD integration gate and config/aging.yml:",
             "region_reading.same_n_members.")
    }
    all_survive <- all(members$survives[members$fitted])

    reading <- if (!primary_supported) {
        if (prim$p < as.numeric(rr$primary_alpha)) "OPPOSITE_DIRECTION" else "NOT_SUPPORTED"
    } else if (all_survive) {
        "SUPPORTED_SURVIVES_GATING_SENSITIVITIES"
    } else {
        "PRIMARY_ONLY_FAILS_GATING_SENSITIVITY"
    }
    list(members = members, prim = prim, primary_supported = primary_supported,
         all_survive = all_survive, reading = reading)
}
