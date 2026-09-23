#### Stable feature identity for the transcription/splicing assays ####
##
## WHY THIS FILE EXISTS
##
## Until 2026-09-23 the PSI (splicing) side of this module keyed features on
## `psi_uid`, the identifier carried by `psi-annotation.tsv`. That identifier is
## POSITIONAL: it is literally `p{row index - 1}` in the order the annotation
## file happens to be written. Verified on 2026-09-23 for all three regions --
## every row satisfies `psi_uid == paste0("p", .I - 1L)`.
##
## The three regions each ship their own `psi-annotation.tsv`. All three hold the
## same 690,907 events (the sorted `event_info` column has one SHA across the
## three files) but in DIFFERENT row orders, so the same `psi_uid` names a
## different splicing event in each region:
##
##     caudate     vs dlpfc        same event for  93,580 / 690,907 psi_uids
##     caudate     vs hippocampus  same event for       0 / 690,907 psi_uids
##     dlpfc       vs hippocampus  same event for 309,180 / 690,907 psi_uids
##
## `config/transcription_splicing.yml` names ONE annotation path,
## `inputs/counts/psi-annotation.tsv`, and that path is a tracked git symlink
## into the CAUDATE delivery. Stage 01 therefore built every region's link table
## from caudate's row order, and stage 02 looked those positional ids up in the
## region's own assay. Because the identifier always exists in every region, the
## `%in%` join always succeeded: no error, no warning, no reduced feature count,
## and a silently different splicing event on 86.5% of DLPFC links and 100% of
## hippocampus links.
##
## THE FIX. A feature is keyed on its own biology -- `event_info` (event type
## plus the exon/junction coordinates) together with `gene_id`, which is unique
## across the annotation (690,907/690,907; `event_info` alone is not, because 308
## events are annotated to two genes). Positional ids survive only as the row
## label the assay object happens to use, they are resolved from the region's own
## metadata at the moment the assay is read, and a key that does not resolve is
## an ERROR rather than a fall-through.
##
## A SECOND DEFECT IN THE ASSAY OBJECTS, found while confirming the first.
## Inside every one of the three `rse-psi.*.RData` objects, `rownames(rse)`
## disagrees with `rowData(rse)`: `rowData` is the region's annotation table in
## FILE order (`rowData$psi_uid` is p0, p1, p2, ... and `rowRanges` agrees with
## it column for column), while the rows themselves are labelled with a
## permutation of the same id set (caudate row 1 is named p63964, DLPFC row 1
## p309180, hippocampus row 1 p402760). Exactly one of the two can describe the
## assay values. A chrY event cannot be quantified in a female donor, so which
## one does was settled from the data, on all three regions:
##
##                       mean NA fraction of chrY-labelled rows
##     labelling         males    females   females - males
##     rownames(rse)     0.67     0.90-0.92     +0.221 to +0.229
##     rowData(rse)      0.60     0.60-0.72     +0.004 to +0.049
##
## `rownames(rse)` is the identifier the assay values follow; `rowData(rse)` is
## the misaligned half. So `rowData` is used here ONLY as a keyed lookup table
## (psi_uid -> event), never positionally against the rows, and that restriction
## is what `psi_feature_table_from_rse()` enforces. Read as a table it is exact:
## `rowData(rse)$psi_uid` and `$event_info` are identical to the region's
## `psi-annotation.tsv` columns, which is what makes it a usable arbiter of which
## region an annotation file belongs to.
##
## The expression side was checked and is clean. `rownames(rse_gene)` equals
## `rowData(rse_gene)$gene_id` equals `gene-annotation.tsv`'s `gene_id`, in
## identical order, in all three regions; the gene annotation file is
## byte-identical across the three deliveries (one MD5); and the identifier is a
## versioned Ensembl gene ID, not a row position. The caudate symlink for
## `gene-annotation.tsv` is therefore harmless, and the expression and ABC
## modalities are unaffected by any of the above.

## Feature-identifier namespaces. Written into every link table so a consumer
## cannot mistake one identifier scheme for another, and so a link table built
## under the withdrawn positional scheme is rejected instead of being reused.
PSI_FEATURE_NAMESPACE  <- "psi_event_key_v1"
GENE_FEATURE_NAMESPACE <- "gene_id_v1"

## `event_info` and `gene_id` contain no "|" anywhere in any of the three
## annotation files (checked 2026-09-23), so the separator cannot make two
## distinct events collide.
PSI_KEY_SEP <- "|"

#' Stable identifier for a splicing event: its own type and coordinates, plus
#' the gene it is annotated to.
psi_event_key <- function(event_info, gene_id) {
    paste(as.character(event_info), as.character(gene_id), sep = PSI_KEY_SEP)
}

#' TRUE for an identifier in the withdrawn positional scheme.
#'
#' Used to refuse a legacy link table outright. Under the positional scheme a
#' `%in%` join succeeds for every row whatever region it came from, which is the
#' property that made the defect silent, so the identifier shape itself has to be
#' treated as a fatal condition rather than something to fall back on.
is_positional_psi_id <- function(x) grepl("^p[0-9]+$", as.character(x))

#' Load the single SummarizedExperiment held in an .RData assay file.
#'
#' Factored out of stage 02 so stage 02 and the smoke checks agree on what
#' "the assay object" means.
load_assay_rse <- function(path) {
    if (!file.exists(path)) stop("Assay file not found: ", path)
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    objs <- mget(ls(env), envir = env)
    hit <- which(vapply(objs, function(x)
        methods::is(x, "SummarizedExperiment"), logical(1)))
    if (length(hit) == 0) {
        stop("No SummarizedExperiment in ", path, "; found: ",
             paste(vapply(objs, function(x) class(x)[1], character(1)),
                   collapse = ", "))
    }
    objs[[hit[1]]]
}

#' The region's own PSI annotation file.
#'
#' `config/transcription_splicing.yml` is `pi_locked` and names a single
#' `annotation.psi` path, which is a symlink into the caudate delivery -- the
#' immediate cause of the defect above. The path is therefore derived from the
#' per-region entry the config DOES carry, `assay_files.psi.{region}`: the
#' annotation that describes an assay is delivered in the same directory as that
#' assay. Nothing is hard-coded here, and there is deliberately no fall-back to
#' `annotation.psi`: reading another region's positional table is the failure
#' being removed, so its absence must stop the run.
psi_annotation_path <- function(ts, region, root = repo_root()) {
    rel <- ts$assay_files$psi[[region]]
    if (is.null(rel) || !nzchar(rel)) {
        stop("config/transcription_splicing.yml has no assay_files.psi entry ",
             "for region '", region, "'")
    }
    assay_abs <- file.path(root, rel)
    if (!file.exists(assay_abs)) {
        stop("PSI assay for ", region, " not found: ", assay_abs)
    }
    ## normalizePath follows the inputs/counts/ symlink into the delivery.
    beside <- file.path(dirname(normalizePath(assay_abs, mustWork = TRUE)),
                        "psi-annotation.tsv")
    if (!file.exists(beside)) {
        stop("No psi-annotation.tsv beside the ", region, " PSI assay.\n",
             "  assay    : ", assay_abs, "\n",
             "  expected : ", beside, "\n",
             "  config annotation.psi (", ts$annotation$psi, ") is NOT used as ",
             "a fall-back: it is a symlink to one region's copy, and psi_uid is ",
             "a row position, so another region's table silently renames every ",
             "event. Deliver the region's annotation beside its assay, or have ",
             "the PI make annotation.psi a per-region map.")
    }
    beside
}

PSI_ANNOT_COLS <- c("psi_uid", "event_info", "event_type", "gene_id",
                    "chrom", "start", "end")

#' Validate a PSI annotation table and attach the stable event key.
#'
#' @param d data.table or data.frame with at least PSI_ANNOT_COLS.
#' @param source_label what to name in an error message.
psi_feature_table <- function(d, source_label = "PSI annotation") {
    d <- data.table::as.data.table(d)
    miss <- setdiff(PSI_ANNOT_COLS, names(d))
    if (length(miss)) {
        stop(source_label, " lacks required column(s): ",
             paste(miss, collapse = ", "))
    }
    out <- d[, .(psi_uid = as.character(psi_uid),
                 event_info = as.character(event_info),
                 event_type = as.character(event_type),
                 gene_id = as.character(gene_id),
                 chrom = as.character(chrom),
                 start = as.integer(start),
                 end = as.integer(end))]
    out[, event_key := psi_event_key(event_info, gene_id)]
    if (anyDuplicated(out$psi_uid)) {
        stop(source_label, ": psi_uid is not unique (",
             sum(duplicated(out$psi_uid)), " duplicates)")
    }
    ## The whole point of the re-key is that this identifier is one-to-one with a
    ## splicing event. If it is not, the join would be ambiguous and we are back
    ## to picking a row by accident.
    if (anyDuplicated(out$event_key)) {
        dup <- out$event_key[duplicated(out$event_key)]
        stop(source_label, ": event_info + gene_id is not unique (",
             length(dup), " collisions), so it cannot key the assay. First: ",
             dup[1])
    }
    data.table::setkeyv(out, "event_key")
    out
}

#' The region's authoritative psi_uid -> event map, taken from its own assay.
#'
#' `rowData(rse)` is read as a KEYED TABLE and never positionally against the
#' rows: in the delivered objects its row order does not match the assay's (see
#' the header). Keying by psi_uid is correct either way -- it stays correct if
#' the objects are ever rebuilt with the two halves aligned.
#'
#' `rownames(rse)` must be a permutation of the same id set, because that is what
#' makes the resolved id usable to select assay rows. A mismatch means the
#' object's labels and metadata describe different feature sets and there is no
#' safe way to proceed.
psi_feature_table_from_rse <- function(rse, label = "PSI assay rowData") {
    map <- psi_feature_table(
        as.data.frame(SummarizedExperiment::rowData(rse)), label)
    rn <- rownames(rse)
    if (is.null(rn)) stop(label, ": the assay has no rownames to resolve")
    if (!setequal(rn, map$psi_uid)) {
        stop(label, ": rownames(rse) is not a permutation of rowData's psi_uid ",
             "(", length(setdiff(rn, map$psi_uid)), " rownames absent from ",
             "rowData, ", length(setdiff(map$psi_uid, rn)), " the other way). ",
             "The object's row labels and its metadata describe different ",
             "feature sets; the assay cannot be keyed.")
    }
    ## Recorded, not fatal: every delivered object has this defect, and the chrY
    ## test above establishes that the ROW LABELS are the half the assay values
    ## follow, so a keyed read is sound. It is surfaced so a run states the
    ## condition of its input instead of depending on it silently.
    attr(map, "rowdata_aligned_with_rownames") <- identical(
        as.character(rn),
        as.character(SummarizedExperiment::rowData(rse)$psi_uid))
    map
}

#' Fail unless the annotation file a link table was built from is this region's.
#'
#' This is the check that makes the defect impossible rather than merely fixed.
#' `rowData` read as a table is exactly the region's own `psi-annotation.tsv`, so
#' the assay object is the arbiter of which region a file belongs to. If stage 01
#' read caudate's copy for a DLPFC run, 86.5% of psi_uids map to a different
#' event key here and the run stops -- including in the case where stage 01 and
#' stage 02 resolve the same wrong path, which a checksum comparison between the
#' two stages would not catch.
psi_assert_annotation_matches_assay <- function(map, annot_path, region) {
    ann <- psi_feature_table(data.table::fread(annot_path),
                             paste0("psi-annotation.tsv for ", region))
    shared <- intersect(ann$psi_uid, map$psi_uid)
    if (length(shared) == 0) {
        stop("The PSI annotation at ", annot_path, " shares no psi_uid with ",
             "the ", region, " assay; it does not describe this assay.")
    }
    a <- ann$event_key[match(shared, ann$psi_uid)]
    b <- map$event_key[match(shared, map$psi_uid)]
    n_bad <- sum(a != b)
    if (n_bad > 0) {
        i <- which(a != b)[1]
        stop("The PSI annotation at ", annot_path, " is NOT the ", region,
             " annotation: ", n_bad, "/", length(shared), " psi_uids name a ",
             "different event in the assay.\n",
             "  psi_uid ", shared[i], " is\n",
             "    annotation: ", a[i], "\n",
             "    assay     : ", b[i], "\n",
             "  psi_uid is a ROW POSITION, and the three regions ship the same ",
             "690,907 events in different orders, so a table from the wrong ",
             "region renames events without failing to join.")
    }
    invisible(list(n_shared = length(shared),
                   n_psi_uid_annotation = nrow(ann),
                   n_psi_uid_assay = nrow(map)))
}

#' Resolve stable event keys to the assay row names that carry them.
#'
#' Every key must resolve. A link table and an assay that disagree about the
#' feature universe is a wrong input, not a subset to be quietly analysed: the
#' silent version of exactly that is the defect this file documents.
psi_resolve_rows <- function(map, keys, context = "link table") {
    keys <- unique(as.character(keys))
    if (length(keys) == 0) stop(context, " declares no PSI feature")
    if (any(is_positional_psi_id(keys))) {
        n <- sum(is_positional_psi_id(keys))
        stop(context, " keys features on the withdrawn positional psi_uid ",
             "scheme (", n, "/", length(keys), " look like p<number>, e.g. ",
             keys[is_positional_psi_id(keys)][1], "). Those ids are row ",
             "positions in one region's annotation file and mean a different ",
             "event in another region. Rebuild the links with stage 01.")
    }
    i <- match(keys, map$event_key)
    if (anyNA(i)) {
        bad <- keys[is.na(i)]
        stop(context, ": ", length(bad), "/", length(keys), " PSI event keys do ",
             "not exist in this region's assay. First unresolved: ", bad[1],
             "\n  The link table and the assay describe different event ",
             "universes (a different annotation version, or a link table built ",
             "for another cell).")
    }
    map$psi_uid[i]
}

#' Resolve gene-modality feature ids, with the same all-or-nothing contract.
#'
#' The gene assay's rownames are its `gene_id`s and equal `gene-annotation.tsv`
#' exactly in all three regions, so a shortfall here is a real annotation/assay
#' version mismatch rather than ordinary attrition.
gene_resolve_rows <- function(rse, keys, context = "link table") {
    keys <- unique(as.character(keys))
    if (length(keys) == 0) stop(context, " declares no gene feature")
    missing <- setdiff(keys, rownames(rse))
    if (length(missing)) {
        stop(context, ": ", length(missing), "/", length(keys), " gene ids are ",
             "absent from the assay. First: ", missing[1],
             "\n  The gene annotation and the gene assay are different ",
             "versions; they must be the same build.")
    }
    keys
}
