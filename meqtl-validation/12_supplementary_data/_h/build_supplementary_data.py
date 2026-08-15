#!/usr/bin/env python3
"""Build sanitized, deterministic Supplementary Data 10-15 archives.

The archives contain derived, manuscript-facing tables only. The builder rejects
duplicate headers, individual-level identifier columns, and cluster-local paths.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import os
from pathlib import Path
import shutil
import tempfile
import zipfile


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT = PROJECT_ROOT / "meqtl-validation/12_supplementary_data/_m"
FIXED_ZIP_TIME = (2026, 8, 7, 0, 0, 0)
REQUIRED_ANALYSIS_SCHEMA_VERSION = 2

REGIONS = ("caudate", "dlpfc", "hippocampus")

PACKAGE_INFO = {
    10: (
        "internal_cpg_cis_meqtl",
        "Internal CpG-level cis-meQTL results and analysis QC",
    ),
    11: (
        "vmr_meqtl_burden",
        "VMR-level aggregation of internal CpG cis-meQTL evidence",
    ),
    12: (
        "external_brain_meqtl_validation",
        "External adult-brain CpG meQTL validation",
    ),
    13: (
        "cross_region_and_donor_group_meqtl",
        "Cross-region and cross-donor-group genetic architecture",
    ),
    14: (
        "transcriptional_and_technical_validation",
        "Transcriptional coupling and technical robustness",
    ),
    15: (
        "schizophrenia_risk_application",
        "Focused schizophrenia-risk regulatory application",
    ),
}

README_TEXT = {
    10: """# Supplementary Data 10: Internal CpG-level cis-meQTL results

This archive provides one permutation-based lead cis-meQTL result per tested CpG,
separately for caudate, DLPFC, and hippocampus and for the Black American (AA in
analysis filenames) and white American (EA in analysis filenames) donor groups.
The primary Black American model is the prespecified M3a covariate specification.
Files include effect estimates, standard errors, allele frequency, CpG-variant
distance, nominal and permutation p-values, and region-specific q-values.

The `audit/` directory contains aggregate sample-overlap, covariate, genome-build,
and external-resource dictionaries. No individual-level methylation, genotype,
covariate, donor identifier, or sample identifier is included.
""",
    11: """# Supplementary Data 11: VMR CpG meQTL burden

This archive aggregates internal CpG-level cis-meQTL evidence to predefined VMRs.
It contains VMR-level burden tables, adjusted burden models, matched high-versus-low
predictability comparisons, and aggregation summaries for three brain regions and
both donor groups. Caudate EA M3a results are retained as a labeled sensitivity
analysis. Duplicate source columns are removed during packaging.

`local_predictability` is an elastic-net cross-validated local SNP-predictability
measure. It is not a calibrated locus-level heritability estimate.
""",
    12: """# Supplementary Data 12: External brain meQTL validation

This archive contains harmonization/overlap summaries, VMR-level external support,
and complete-assay-universe adjusted/matched analyses for Jaffe DLPFC. Schulz
hippocampus is retained as positive-only supporting overlap because tested negatives
are unavailable. BrainSeq/LIBD results are explicitly exploratory because the cohort
overlaps the internal discovery resource and is not independent validation.

Original third-party supplementary tables and liftOver chain files are not
redistributed. Consult `audit/public_meqtl_resources.tsv` in Supplementary Data 10
for citations, access locations, genome builds, and inclusion roles.
""",
    13: """# Supplementary Data 13: Cross-region and donor-group architecture

This archive contains cross-region CpG effect concordance, VMR sharing, shared-donor
genotype-by-region results, caudate sample-size downsampling summaries, Black
American-white American CpG effect concordance, predictability portability, and
MAF/local-SNP-opportunity matched comparisons.

Donor and sample lists are excluded. Donor-group contrasts are interpreted as
portability and genetic-context comparisons, not as ancestry-specific effects.
""",
    14: """# Supplementary Data 14: Transcriptional and technical validation

This archive links VMR CpG meQTL support to existing expression and splicing
associations and provides repeat/chromatin robustness analyses. It includes the
consolidated LINE/L1, H3K9me3, quiescent-chromatin, mappability, SNP-proximity,
segmental-duplication, and cell-composition sensitivity results.

Only VMR-level aggregate cell-composition metrics are included. Sample-level cell
covariates and sample-by-VMR methylation matrices are excluded.
""",
    15: """# Supplementary Data 15: Schizophrenia-risk application

This archive supports the focused schizophrenia-risk regulatory application. It
contains prespecified risk-locus definitions, risk-variant-CpG tests, VMR links,
architecture models, regional comparisons, downsampling, transcriptional coupling,
targeted external eQTL support, diagnosis sensitivity analyses, and tables underlying
the prioritized locus panels.

The results support regulatory association, not mediation, colocalization, or
causality. Raw PGC/GTEx/GWAS source files, internal retain/omit decisions, claim
snapshots, figures, donor lists, and regional GWAS slices are not redistributed.
""",
}

SENSITIVE_COLUMNS = {
    "sample_id",
    "sampleid",
    "donor_id",
    "donorid",
    "subject_id",
    "individual_id",
    "brnum",
    "rnum",
}
LOCAL_PATH_MARKERS = (b"/projects/", b"/gpfs/", b"/home/", b"/scratch/")
KNOWN_PATH_REPLACEMENTS = {
    "/projects/b1213/users/kynon/projects/dna-methylation-heritability/": "",
    "/gpfs/projects/b1213/users/kynon/projects/dna-methylation-heritability/": "",
    "/projects/b1213/users/alexis/projects/dna-methylation-heritability/": "external_staff_repository/",
    "/gpfs/projects/b1213/users/alexis/projects/dna-methylation-heritability/": "external_staff_repository/",
    "/projects/b1213/resources/": "institutional_resource/",
    "/gpfs/projects/b1213/resources/": "institutional_resource/",
}


def add(files: list[tuple[str, str]], source: str, target: str | None = None) -> None:
    files.append((source, target or source))


def require_repaired_outputs() -> None:
    """Fail closed rather than package pre-repair manuscript results."""
    failures: list[str] = []
    for region in REGIONS:
        path = (
            PROJECT_ROOT
            / "meqtl-validation/02_vmr_meqtl_burden/_m"
            / region
            / "vmr_meqtl_burden.tsv.gz"
        )
        if not path.exists():
            failures.append(f"missing {path}")
            continue
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            first = next(reader, None)
        try:
            version = int(float((first or {}).get("analysis_schema_version", "")))
        except ValueError:
            version = -1
        if version < REQUIRED_ANALYSIS_SCHEMA_VERSION:
            failures.append(
                f"{path} has schema {version}; need {REQUIRED_ANALYSIS_SCHEMA_VERSION}"
            )
    if failures:
        raise SystemExit(
            "Refusing to package stale meQTL-validation outputs:\n- "
            + "\n- ".join(failures)
            + "\nRerun in the order documented in meqtl-validation/REPAIR_V2.md."
        )


def package_files(number: int) -> list[tuple[str, str]]:
    files: list[tuple[str, str]] = []

    if number == 10:
        for name in (
            "sample_overlap.tsv",
            "covariate_dictionary.tsv",
            "genome_build_audit.tsv",
            "public_meqtl_resources.tsv",
        ):
            add(files, f"inputs/data_dictionary/_m/{name}", f"audit/{name}")
        for region in REGIONS:
            base = f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m"
            add(files, f"{base}/preflight/preflight_summary.tsv", f"AA/{region}/preflight_summary.tsv")
            add(files, f"{base}/covariate_sensitivity/comparison_summary.tsv", f"AA/{region}/covariate_model_comparison.tsv")
            add(files, f"{base}/covariate_sensitivity/covariate_model_manifest.tsv", f"AA/{region}/covariate_model_manifest.tsv")
            add(files, f"{base}/tensorqtl/cpg_meqtl_{region}.cis_qtl.txt.gz", f"AA/{region}/lead_cis_meqtl_per_cpg.tsv.gz")
            add(files, f"{base}/tensorqtl/qc/meqtl_qc_summary.tsv", f"AA/{region}/meqtl_qc_summary.tsv")
            add(files, f"{base}/tensorqtl/qc/pvalue_histogram.tsv", f"AA/{region}/pvalue_histogram.tsv")
            add(files, f"{base}/tensorqtl/qc/calibration/calibration_summary.tsv", f"AA/{region}/calibration_summary.tsv")
            add(files, f"{base}/tensorqtl/qc/calibration/lambda_by_qval.tsv", f"AA/{region}/lambda_by_qval.tsv")
            add(files, f"{base}/preflight/EA/preflight_summary.tsv", f"EA/{region}/preflight_summary.tsv")
            add(files, f"{base}/tensorqtl/EA/cpg_meqtl_{region}_EA.cis_qtl.txt.gz", f"EA/{region}/lead_cis_meqtl_per_cpg.tsv.gz")
            add(files, f"{base}/tensorqtl/EA/qc/meqtl_qc_summary.tsv", f"EA/{region}/meqtl_qc_summary.tsv")
            add(files, f"{base}/tensorqtl/EA/qc/pvalue_histogram.tsv", f"EA/{region}/pvalue_histogram.tsv")

    elif number == 11:
        for group, source_group in (("AA", ""), ("EA", "EA/")):
            for region in REGIONS:
                base = f"meqtl-validation/02_vmr_meqtl_burden/_m/{source_group}{region}"
                for name in (
                    "vmr_meqtl_burden.tsv.gz",
                    "aggregation_summary.tsv",
                    "burden_model_results.tsv",
                    "matched_analysis_results.tsv",
                ):
                    add(files, f"{base}/{name}", f"{group}/{region}/{name}")
        base = "meqtl-validation/02_vmr_meqtl_burden/_m/EA_M3a/caudate"
        for name in (
            "vmr_meqtl_burden.tsv.gz",
            "aggregation_summary.tsv",
            "burden_model_results.tsv",
            "matched_analysis_results.tsv",
        ):
            add(files, f"{base}/{name}", f"sensitivity/EA_M3a/caudate/{name}")

    elif number == 12:
        base = "meqtl-validation/03_external_meqtl_validation/_m"
        for name in (
            "external_support_models_all.tsv",
            "external_matched_all.tsv",
            "external_support_primary_summary.tsv",
        ):
            add(files, f"{base}/{name}", f"summary/{name}")
        resources = (
            "brainseq_wgbs_meqtl_scz_subset",
            "jaffe_dlpfc_450k_meqtl",
            "schulz_hippocampus_array_meqtl",
        )
        for resource in resources:
            add(files, f"{base}/harmonized/{resource}.overlap_summary.tsv", f"harmonization/{resource}.overlap_summary.tsv")
        for region in REGIONS:
            for resource in resources:
                add(files, f"{base}/{region}/vmr_external_support_{resource}.tsv.gz", f"{region}/{resource}/vmr_external_support.tsv.gz")
                add(files, f"{base}/{region}/external_support_model_{resource}.tsv", f"{region}/{resource}/model_results.tsv")
                add(files, f"{base}/{region}/external_matched_{resource}.tsv", f"{region}/{resource}/matched_results.tsv")

    elif number == 13:
        cross = "meqtl-validation/04_cross_region_sharing/_m"
        for name in (
            "burden_gradient_by_region.tsv",
            "cross_region_cpg_concordance.tsv",
            "cross_region_vmr_sharing.tsv",
            "discovery_power_summary.tsv",
            "vmr_cross_region_support.tsv.gz",
        ):
            add(files, f"{cross}/{name}", f"cross_region/{name}")
        for name in (
            "gxregion_by_predictability_tercile.tsv",
            "gxregion_pair_results.tsv.gz",
            "gxregion_summary.tsv",
            "pairs_tested_definition.tsv.gz",
        ):
            add(files, f"{cross}/gxregion/{name}", f"cross_region/gxregion/{name}")
        down = "meqtl-validation/10_downsampling_caudate/_m"
        for name in (
            "downsample_design_summary.tsv",
            "tensorqtl_downsample_replicate_results.tsv",
            "tensorqtl_downsample_vs_regions.tsv",
        ):
            add(files, f"{down}/{name}", f"cross_region/caudate_downsample/{name}")
        donor = "meqtl-validation/05_donor_group_comparison/_m"
        for name in (
            "aa_ea_effect_concordance_summary.tsv",
            "aa_ea_predictability_portability.tsv",
            "concordant_high_genomic_annotation.tsv",
            "donor_group_coefficient_comparison.tsv",
            "maf_ld_discovery_strata.tsv",
            "maf_ld_matched_discovery.tsv",
        ):
            add(files, f"{donor}/{name}", f"donor_group/{name}")
        for region in REGIONS:
            for name in (
                "aa_ea_cpg_effect_concordance.tsv.gz",
                "aa_ea_predictability_joined.tsv.gz",
                "concordant_high_predictability_vmrs.tsv.gz",
                "donor_group_coefficient_comparison.tsv",
            ):
                add(files, f"{donor}/{region}/{name}", f"donor_group/{region}/{name}")

    elif number == 14:
        tx = "meqtl-validation/06_transcription_splicing_integration/_m"
        for name in (
            "tx_enrichment_all.tsv",
            "tx_enrichment_primary.tsv",
            "tx_enrichment_both_modalities.tsv",
        ):
            add(files, f"{tx}/{name}", f"transcription/{name}")
        for region in REGIONS:
            for name in (
                "meqtl_x_expression_enrichment.tsv",
                "meqtl_x_expression_abc_enrichment.tsv",
                "meqtl_x_psi_enrichment.tsv",
                "vmr_meqtl_expression_joined.tsv.gz",
                "vmr_meqtl_expression_abc_joined.tsv.gz",
                "vmr_meqtl_psi_joined.tsv.gz",
            ):
                add(files, f"{tx}/{region}/{name}", f"transcription/{region}/{name}")
        tech = "meqtl-validation/07_repeat_mappability_sensitivity/_m"
        for name in (
            "annotation_asset_manifest.tsv",
            "consolidated_robustness_table.tsv",
            "tech_join_completeness.tsv",
            "vmr_technical_annotation_index.tsv",
        ):
            add(files, f"{tech}/{name}", f"technical_robustness/{name}")
        for region in REGIONS:
            for name in (
                "burden_tech_join.tsv.gz",
                "robustness_results.tsv",
                "vmr_features_with_tech.tsv.gz",
                "vmr_technical_annotations.tsv",
            ):
                add(files, f"{tech}/{region}/{name}", f"technical_robustness/{region}/{name}")
        cell = "meqtl-validation/11_celltype_compartment_sensitivity/_m"
        for name in (
            "burden_celltype_models.tsv",
            "celltype_robustness_table.tsv",
            "celltype_sample_overlap.tsv",
            "consolidated_robustness_table_with_celltype.tsv",
            "enrichment_celltype_models.tsv",
            "sample_level_celltype_checks.tsv",
            "vmr_cell_metrics_summary.tsv",
        ):
            add(files, f"{cell}/{name}", f"cell_composition/{name}")
        for region in REGIONS:
            for name in ("vmr_cell_metrics.tsv.gz", "vmr_features_with_cell_metrics.tsv.gz"):
                add(files, f"{cell}/{region}/{name}", f"cell_composition/{region}/{name}")

    elif number == 15:
        scz = "meqtl-validation/08_schizophrenia_risk_application/_m"
        for region in REGIONS:
            for name in (
                "architecture_matched_permutation.tsv",
                "architecture_model_results.tsv",
                "define_loci_summary.tsv",
                "link_summary.tsv",
                "risk_variant_cpg_failures.tsv",
                "risk_variant_cpg_meqtl.tsv.gz",
                "risk_variant_cpg_tests_summary.tsv",
                "scz_index_snps_hg38.tsv",
                "scz_meqtl_tx_coupled_vmrs.tsv.gz",
                "scz_risk_loci_hg38.tsv",
                "scz_vmr_level_summary.tsv.gz",
                "tx_integration_results.tsv",
                "vmr_locus_links.tsv.gz",
            ):
                add(files, f"{scz}/{region}/{name}", f"region_results/{region}/{name}")
        for name in (
            "discovery_power_summary.tsv",
            "locus_by_region.tsv.gz",
            "pairwise_concordance.tsv",
            "sig_pairs_by_region.tsv.gz",
        ):
            add(files, f"{scz}/cross_region/{name}", f"cross_region/{name}")
        for name in (
            "gxregion_pair_results.tsv.gz",
            "gxregion_prioritized_locus_summary.tsv",
            "gxregion_summary.tsv",
            "pairs_tested_definition.tsv.gz",
        ):
            add(files, f"{scz}/gxregion/{name}", f"cross_region/gxregion/{name}")
        for name in (
            "downsample_design_summary.tsv",
            "downsample_prioritized_locus_results.tsv.gz",
            "downsample_prioritized_stability.tsv",
            "downsample_replicate_results.tsv",
            "downsample_run_summary.tsv",
            "downsample_vs_regions.tsv",
        ):
            add(files, f"{scz}/caudate_downsample/{name}", f"downsampling/{name}")
        for name in (
            "diagnosis_association_results.tsv",
            "diagnosis_locus_rollup.tsv",
            "diagnosis_sensitivity_stability.tsv",
            "tested_transcript_features.tsv",
            "tested_vmrs.tsv",
            "vmr_mean_methylation_meta.tsv",
        ):
            add(files, f"{scz}/diagnosis/{name}", f"diagnosis/{name}")
        for name in (
            "level3_gene_targets.tsv",
            "level3_gene_targets_summary.tsv",
            "level3_locus_summary.tsv",
            "level3_unique_genes.tsv",
            "level3_vmr_feature_links.tsv",
            "libd_risk_variant_eqtl.tsv.gz",
            "libd_risk_variant_eqtl_summary.tsv",
        ):
            add(files, f"{scz}/level3/{name}", f"shared_genetic_support/{name}")
        for name in ("gtex_level3_locus_gene_summary.tsv", "gtex_level3_signif_hits.tsv.gz"):
            add(files, f"{scz}/level3/gtex/{name}", f"shared_genetic_support/gtex/{name}")
        for name in (
            "locus_priority_ranked_all.tsv.gz",
            "prioritization_summary.tsv",
            "prioritized_loci.tsv",
            "prioritized_loci_rationale.tsv",
        ):
            add(files, f"{scz}/prioritized/{name}", f"prioritized_loci/{name}")
        add(files, f"{scz}/locus_panels/locus_panel_manifest.tsv", "prioritized_loci/locus_panel_manifest.tsv")
        panel_root = PROJECT_ROOT / scz / "locus_panels"
        for locus_dir in sorted(path for path in panel_root.iterdir() if path.is_dir() and path.name != "figures"):
            for name in (
                "cross_region_forest.tsv.gz",
                "gtex_level3_hits.tsv",
                "locus_meta.tsv",
                "meqtl_caudate_cpgs.tsv.gz",
                "tx_links_fdr.tsv",
                "vmr_predictability.tsv",
            ):
                path = locus_dir / name
                if path.exists():
                    add(files, str(path.relative_to(PROJECT_ROOT)), f"prioritized_loci/panel_data/{locus_dir.name}/{name}")
    else:
        raise ValueError(f"Unknown package number: {number}")

    return files


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def table_header(path: Path) -> list[str]:
    if not (path.name.endswith(".tsv") or path.name.endswith(".tsv.gz") or path.name.endswith(".txt.gz")):
        return []
    with open_text(path) as handle:
        line = handle.readline().rstrip("\r\n")
    return line.split("\t") if line else []


def remove_duplicate_columns(source: Path, target: Path) -> None:
    """Copy a TSV while retaining only the first occurrence of each header."""
    target.parent.mkdir(parents=True, exist_ok=True)
    with open_text(source) as src:
        header = src.readline().rstrip("\r\n").split("\t")
        keep: list[int] = []
        seen: set[str] = set()
        for index, column in enumerate(header):
            if column not in seen:
                seen.add(column)
                keep.append(index)

        if target.suffix == ".gz":
            raw = target.open("wb")
            gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
            out = io.TextIOWrapper(gz, encoding="utf-8", newline="")
        else:
            raw = None
            gz = None
            out = target.open("wt", encoding="utf-8", newline="")
        try:
            out.write("\t".join(header[i] for i in keep) + "\n")
            for line in src:
                fields = line.rstrip("\r\n").split("\t")
                out.write("\t".join(fields[i] if i < len(fields) else "" for i in keep) + "\n")
        finally:
            out.close()
            if gz is not None:
                gz.close()
            if raw is not None:
                raw.close()


def copy_source(source: Path, target: Path) -> None:
    if not source.exists():
        raise FileNotFoundError(source)
    header = table_header(source)
    if header and len(header) != len(set(header)):
        remove_duplicate_columns(source, target)
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    if target.suffix != ".gz" and table_header(target):
        text = target.read_text(encoding="utf-8")
        sanitized = text
        for old, new in KNOWN_PATH_REPLACEMENTS.items():
            sanitized = sanitized.replace(old, new)
        if sanitized != text:
            target.write_text(sanitized, encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_rows(path: Path) -> int | str:
    if not table_header(path):
        return "NA"
    with open_text(path) as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def validate_table(path: Path) -> None:
    header = table_header(path)
    if not header:
        return
    duplicates = sorted({name for name in header if header.count(name) > 1})
    if duplicates:
        raise ValueError(f"Duplicate columns in {path}: {duplicates}")
    normalized = {name.lower().replace("-", "_") for name in header}
    sensitive = sorted(normalized & SENSITIVE_COLUMNS)
    if sensitive:
        raise ValueError(f"Individual-level identifier columns in {path}: {sensitive}")

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            encoded = line.encode("utf-8", errors="ignore")
            if any(marker in encoded for marker in LOCAL_PATH_MARKERS):
                raise ValueError(f"Cluster-local path in {path} line {line_number}")


def write_schema(package_root: Path, data_files: list[Path]) -> None:
    schema = package_root / "SCHEMA.tsv"
    with schema.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("file", "column_index", "column_name"))
        for path in sorted(data_files):
            for index, name in enumerate(table_header(path), start=1):
                writer.writerow((str(path.relative_to(package_root)), index, name))


def write_manifest(package_root: Path) -> None:
    manifest = package_root / "file_manifest.tsv"
    paths = sorted(path for path in package_root.rglob("*") if path.is_file() and path != manifest)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("file", "bytes", "sha256", "n_data_rows"))
        for path in paths:
            writer.writerow((str(path.relative_to(package_root)), path.stat().st_size, sha256(path), count_rows(path)))


def make_zip(package_root: Path, output_zip: Path) -> None:
    output_zip.parent.mkdir(parents=True, exist_ok=True)
    temp_zip = output_zip.with_suffix(".zip.tmp")
    with zipfile.ZipFile(temp_zip, "w") as archive:
        for path in sorted(package_root.rglob("*")):
            if not path.is_file():
                continue
            arcname = f"{package_root.name}/{path.relative_to(package_root)}"
            info = zipfile.ZipInfo(arcname, FIXED_ZIP_TIME)
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            compression = zipfile.ZIP_STORED if path.suffix == ".gz" else zipfile.ZIP_DEFLATED
            archive.writestr(info, path.read_bytes(), compress_type=compression, compresslevel=9)
    os.replace(temp_zip, output_zip)


def build_package(number: int, staging_root: Path, output_root: Path) -> tuple[Path, int]:
    slug, _description = PACKAGE_INFO[number]
    package_root = staging_root / f"supplementary_data_{number}_{slug}"
    package_root.mkdir(parents=True)
    data_files: list[Path] = []
    for source_rel, target_rel in package_files(number):
        target = package_root / target_rel
        copy_source(PROJECT_ROOT / source_rel, target)
        validate_table(target)
        data_files.append(target)

    (package_root / "README.md").write_text(README_TEXT[number].strip() + "\n", encoding="utf-8")
    write_schema(package_root, data_files)
    write_manifest(package_root)

    output_zip = output_root / package_root.name / f"{slug}.zip"
    make_zip(package_root, output_zip)
    return output_zip, len(data_files)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require_repaired_outputs()
    output_root = args.output.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[object, ...]] = []
    with tempfile.TemporaryDirectory(prefix="meqtl_sdata_", dir="/tmp") as temp:
        staging_root = Path(temp)
        for number in sorted(PACKAGE_INFO):
            output_zip, n_files = build_package(number, staging_root, output_root)
            _slug, description = PACKAGE_INFO[number]
            rows.append((number, output_zip.relative_to(output_root), output_zip.stat().st_size, sha256(output_zip), n_files, description))

    manifest = output_root / "supplementary_data_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("supplementary_data", "archive", "bytes", "sha256", "n_data_files", "description"))
        writer.writerows(rows)


if __name__ == "__main__":
    main()
