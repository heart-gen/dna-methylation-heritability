#!/usr/bin/env python3
"""Generate Phase 0 data-dictionary tables for AJHG strengthening analyses.

Extends the meqtl-validation data audit into the formal deliverables required by
AGENTS.md / analysis_strategy.md. Does not fit meQTL or phenotype models.
"""

from __future__ import annotations

import csv
import glob
import os, sys
import platform
from pathlib import Path
from datatime import datetime, timezone
from collections import Counter, defaultdict

try:
    import yaml
except ImportError:  # pragma: no cover - fallback without PyYAML
    yaml = None


SCRIPT_DIR = Path(__file__).resolve().parent
OUTDIR = SCRIPT_DIR.parent / "_m"
PROJECT_ROOT = SCRIPT_DIR.resolve().parents[2]
PATHS_YML = PROJECT_ROOT / "config" / "paths.yml"


def load_paths() -> dict:
    if yaml is not None and PATHS_YML.exists():
        with PATHS_YML.open() as handle:
            return yaml.safe_load(handle)
    # Minimal fallback if PyYAML unavailable
    return {
        "project_root": str(PROJECT_ROOT),
        "staff_project_root": "/projects/b1213/users/alexis/projects/dna-methylation-heritability",
        "regions": ["caudate", "dlpfc", "hippocampus"],
        "phenotype_table": "sample_summary/_m/phenotype_data.tsv",
        "vmr_bed_template": "vmr-analysis/{region}/_m/vmr.bed",
        "local_predictability_summary_template": (
            "heritability/elastic_net_model/all_individuals/{region}/_m/"
            "{region}_summary_elastic-net_AA.tsv"
        ),
        "cell_proportions_template": "inputs/cell_proportions/_m/music-proportions-{region}.tsv",
        "genotype": {
            "AA": {
                "pgen": "inputs/genotypes/TOPMed_LIBD.AA.pgen",
                "psam": "inputs/genotypes/TOPMed_LIBD.AA.psam",
                "pvar": "inputs/genotypes/TOPMed_LIBD.AA.pvar",
                "eigenvec": "inputs/genotypes/TOPMed_LIBD.AA.eigenvec",
                "eigenval": "inputs/genotypes/TOPMed_LIBD.AA.eigenval",
            }
        },
        "cpg_methylation_root": "/projects/b1213/users/alexis/projects/dna-methylation-heritability",
        "cpg_methylation_matrix_template": "vmr-analysis/{region}/_m/cpg/chr_{chrom}/cpg_meth.phen",
        "regulatory_context_root": (
            "heritability/elastic_net_model/all_individuals/tissue_comparison/"
            "regulatory_context/_m"
        ),
        "support_files": {
            "repeat_masker": "inputs/supportfiles/_m/repeat-masker-hg38.gz",
            "hg38_to_hg19_chain": "inputs/supportfiles/_m/hg38ToHg19.over.chain",
        },
        "vmr_annotation_template": (
            "heritability/elastic_net_model/all_individuals/tissue_comparison/"
            "annotation/_m/{region}_vmr_annotations_hg38.tsv"
        ),
    }


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_table(path: Path, delimiter: str | None = None) -> list[dict]:
    if delimiter is None:
        delimiter = "," if path.suffix == ".csv" else "\t"
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def norm_region(value: str) -> str:
    lower = (value or "").strip().lower()
    if lower == "caudate nucleus":
        return "caudate"
    return lower


def nonmissing(row: dict, col: str) -> bool:
    val = row.get(col, "")
    return val not in ("", "NA", "NaN", "nan", "NULL", "None")


def count_lines(path: Path) -> int | str:
    try:
        with path.open("rb") as handle:
            return sum(1 for _ in handle)
    except OSError as exc:
        return f"ERROR: {exc}"


def header_n_columns(path: Path) -> int | str:
    try:
        with path.open("rb") as handle:
            header = handle.readline()
        if not header:
            return 0
        return header.count(b"\t") + 1
    except OSError as exc:
        return f"ERROR: {exc}"


def write_sample_overlap(cfg: dict) -> None:
    project = Path(cfg["project_root"])
    phen_path = project / cfg["phenotype_table"]
    rows = read_table(phen_path)
    for row in rows:
        row["_region"] = norm_region(row.get("region", ""))

    # Genotype sample IDs (psam may lack a header; column 0 is BrNum / FID)
    psam = project / cfg["genotype"]["AA"]["psam"]
    geno_ids: set[str] = set()
    if psam.exists():
        with psam.open() as handle:
            first = handle.readline().rstrip("\n").split("\t")
            headerish = [h.lstrip("#") for h in first]
            has_header = "FID" in headerish or "IID" in headerish
            if has_header:
                fid_idx = headerish.index("FID") if "FID" in headerish else 0
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) > fid_idx:
                        geno_ids.add(parts[fid_idx])
            else:
                geno_ids.add(first[0])
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if parts:
                        geno_ids.add(parts[0])

    # Cell composition sample IDs by region
    cell_ids: dict[str, set[str]] = {}
    for region in cfg["regions"]:
        cell_path = project / cfg["cell_proportions_template"].format(region=region)
        ids: set[str] = set()
        if cell_path.exists():
            for crow in read_table(cell_path):
                sid = crow.get("sample_id") or crow.get("BrNum") or crow.get("brnum") or ""
                if sid:
                    ids.add(sid)
        cell_ids[region] = ids

    # CpG matrix sample IDs: PLINK phen FID column (BrNum); prefer small chr_22
    staff = Path(cfg["cpg_methylation_root"])
    cpg_ids: dict[str, set[str]] = {}
    for region in cfg["regions"]:
        ids: set[str] = set()
        example = None
        for chrom in ("22", "21", "1"):
            candidate = staff / "vmr-analysis" / region / "_m" / "cpg" / f"chr_{chrom}" / "cpg_meth.phen"
            if candidate.exists():
                example = candidate
                break
        if example is not None:
            with example.open() as handle:
                handle.readline()  # site header
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if parts:
                        ids.add(parts[0])  # FID = BrNum
        cpg_ids[region] = ids

    out_rows = []
    for race in sorted({r.get("race", "") for r in rows}):
        for region in cfg["regions"]:
            subset = [r for r in rows if r.get("race") == race and r["_region"] == region]
            if not subset:
                continue
            brnums = {r.get("brnum", "") for r in subset if r.get("brnum")}
            n_pheno = len(subset)
            n_geno = sum(1 for b in brnums if b in geno_ids) if geno_ids else 0
            n_cell = sum(1 for b in brnums if b in cell_ids.get(region, set()))
            n_cpg = sum(1 for b in brnums if b in cpg_ids.get(region, set()))
            primary_covars = ["agedeath", "sex", "primarydx", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5"]
            n_complete = sum(1 for r in subset if all(nonmissing(r, c) for c in primary_covars))
            dx = Counter(r.get("primarydx", "") for r in subset)
            out_rows.append({
                "race": race,
                "region": region,
                "n_phenotype_samples": n_pheno,
                "n_unique_donors": len(brnums),
                "n_with_genotype_id_match": n_geno,
                "n_with_cell_proportion": n_cell,
                "n_with_cpg_matrix_id_match": n_cpg,
                "n_complete_primary_covariates": n_complete,
                "n_control": dx.get("Control", 0),
                "n_schizophrenia": dx.get("Schizo", 0),
                "has_vmr_bed": (project / cfg["vmr_bed_template"].format(region=region)).exists(),
                "has_predictability_summary": (
                    project / cfg["local_predictability_summary_template"].format(region=region)
                ).exists(),
                "has_cpg_matrices_staff": bool(cpg_ids.get(region)),
            })

    # repeated donor summary
    donor_regions = defaultdict(set)
    for r in rows:
        if r.get("brnum"):
            donor_regions[r["brnum"]].add(r["_region"])
    for n_reg, n_donors in sorted(Counter(len(v) for v in donor_regions.values()).items()):
        out_rows.append({
            "race": "ALL",
            "region": f"repeated_donors_{n_reg}_regions",
            "n_phenotype_samples": "",
            "n_unique_donors": n_donors,
            "n_with_genotype_id_match": "",
            "n_with_cell_proportion": "",
            "n_with_cpg_matrix_id_match": "",
            "n_complete_primary_covariates": "",
            "n_control": "",
            "n_schizophrenia": "",
            "has_vmr_bed": "",
            "has_predictability_summary": "",
            "has_cpg_matrices_staff": "",
        })

    write_tsv(
        OUTDIR / "sample_overlap.tsv",
        out_rows,
        [
            "race", "region", "n_phenotype_samples", "n_unique_donors",
            "n_with_genotype_id_match", "n_with_cell_proportion",
            "n_with_cpg_matrix_id_match", "n_complete_primary_covariates",
            "n_control", "n_schizophrenia", "has_vmr_bed",
            "has_predictability_summary", "has_cpg_matrices_staff",
        ],
    )


def write_covariate_dictionary(cfg: dict) -> None:
    rows = [
        {"variable": "agedeath", "role": "primary_covariate", "description": "Age at death", "source_table": cfg["phenotype_table"], "missingness_policy": "required for primary model"},
        {"variable": "sex", "role": "primary_covariate", "description": "Biological sex", "source_table": cfg["phenotype_table"], "missingness_policy": "required for primary model"},
        {"variable": "primarydx", "role": "primary_covariate", "description": "Primary diagnosis (Control / Schizo)", "source_table": cfg["phenotype_table"], "missingness_policy": "required for primary model"},
        {"variable": "snpPC1", "role": "primary_covariate", "description": "Ancestry PC1", "source_table": cfg["phenotype_table"], "missingness_policy": "required for primary model"},
        {"variable": "snpPC2", "role": "primary_covariate", "description": "Ancestry PC2", "source_table": cfg["phenotype_table"], "missingness_policy": "required for primary model"},
        {"variable": "snpPC3", "role": "primary_covariate", "description": "Ancestry PC3", "source_table": cfg["phenotype_table"], "missingness_policy": "required for primary model"},
        {"variable": "snpPC4", "role": "primary_covariate", "description": "Ancestry PC4", "source_table": cfg["phenotype_table"], "missingness_policy": "required for primary model"},
        {"variable": "snpPC5", "role": "primary_covariate", "description": "Ancestry PC5", "source_table": cfg["phenotype_table"], "missingness_policy": "required for primary model"},
        {"variable": "pmi", "role": "technical_sensitivity", "description": "Postmortem interval", "source_table": cfg["phenotype_table"], "missingness_policy": "sensitivity subset"},
        {"variable": "ph", "role": "technical_sensitivity", "description": "Brain pH", "source_table": cfg["phenotype_table"], "missingness_policy": "sensitivity subset"},
        {"variable": "smoking", "role": "toxicology_sensitivity", "description": "Smoking exposure", "source_table": cfg["phenotype_table"], "missingness_policy": "sensitivity subset"},
        {"variable": "nicotine", "role": "toxicology_sensitivity", "description": "Nicotine", "source_table": cfg["phenotype_table"], "missingness_policy": "sensitivity subset"},
        {"variable": "cocaine", "role": "toxicology_sensitivity", "description": "Cocaine", "source_table": cfg["phenotype_table"], "missingness_policy": "sensitivity subset"},
        {"variable": "ethanol", "role": "toxicology_sensitivity", "description": "Ethanol", "source_table": cfg["phenotype_table"], "missingness_policy": "sensitivity subset"},
        {"variable": "amphetamines", "role": "toxicology_sensitivity", "description": "Amphetamines", "source_table": cfg["phenotype_table"], "missingness_policy": "sensitivity subset"},
        {"variable": "antipsychotics", "role": "clinical_sensitivity", "description": "Antipsychotic exposure", "source_table": cfg["phenotype_table"], "missingness_policy": "sensitivity subset"},
        {"variable": "lifetime_antipsych", "role": "clinical_sensitivity", "description": "Lifetime antipsychotic exposure", "source_table": cfg["phenotype_table"], "missingness_policy": "sensitivity subset"},
        {"variable": "cell_proportions", "role": "sensitivity_covariate", "description": "MuSiC cell-type proportions", "source_table": "inputs/cell_proportions/_m/music-proportions-{region}.tsv", "missingness_policy": "sensitivity; not in primary meQTL model by default"},
        {"variable": "latent_factors", "role": "optional_meqtl_covariate", "description": "Region-specific methylation latent factors", "source_table": "derived during Phase 1 calibration", "missingness_policy": "include only if calibration improves without removing genetic signal"},
        {"variable": "batch", "role": "technical_covariate", "description": "Sequencing or processing batch if available", "source_table": "region phenotype tables / WGBS metadata", "missingness_policy": "include when available without large sample loss"},
    ]
    write_tsv(
        OUTDIR / "covariate_dictionary.tsv",
        rows,
        ["variable", "role", "description", "source_table", "missingness_policy"],
    )


def write_genome_build_audit(cfg: dict) -> None:
    project = Path(cfg["project_root"])
    rows = [
        {
            "dataset": "VMR BED intervals",
            "path_template": cfg["vmr_bed_template"],
            "genome_build": "hg38",
            "coordinate_system": "0-based BED",
            "evidence": "filenames and annotation companions use hg38",
            "liftOver_needed": "no",
            "notes": "Primary analysis build",
        },
        {
            "dataset": "VMR genomic annotations",
            "path_template": cfg.get("vmr_annotation_template", ""),
            "genome_build": "hg38",
            "coordinate_system": "interval annotations",
            "evidence": "files named *_vmr_annotations_hg38*.tsv",
            "liftOver_needed": "no",
            "notes": "",
        },
        {
            "dataset": "Local genetic predictability summaries",
            "path_template": cfg["local_predictability_summary_template"],
            "genome_build": "hg38",
            "coordinate_system": "VMR-level (inherits VMR BED)",
            "evidence": "built from hg38 VMRs",
            "liftOver_needed": "no",
            "notes": "",
        },
        {
            "dataset": "AA TOPMed genotypes",
            "path_template": cfg["genotype"]["AA"]["pvar"],
            "genome_build": "hg38",
            "coordinate_system": "variant positions in pvar",
            "evidence": "TOPMed-imputed LIBD genotypes used throughout manuscript",
            "liftOver_needed": "no",
            "notes": "Confirm chrom naming (chr1 vs 1) during Phase 1 preflight",
        },
        {
            "dataset": "CpG methylation matrices (staff repo)",
            "path_template": cfg["cpg_methylation_matrix_template"],
            "genome_build": "hg38",
            "coordinate_system": "CpG site columns / BED conversion",
            "evidence": "WGBS pipeline aligned to hg38",
            "liftOver_needed": "no",
            "notes": "Matrices live under Alexis repo; project cpg dirs empty",
        },
        {
            "dataset": "RepeatMasker",
            "path_template": cfg.get("support_files", {}).get("repeat_masker", ""),
            "genome_build": "hg38",
            "coordinate_system": "interval",
            "evidence": "filename repeat-masker-hg38.gz",
            "liftOver_needed": "no",
            "notes": "",
        },
        {
            "dataset": "hg38ToHg19 chain",
            "path_template": cfg.get("support_files", {}).get("hg38_to_hg19_chain", ""),
            "genome_build": "chain",
            "coordinate_system": "liftOver chain",
            "evidence": "file present for S-LDSC / legacy analyses",
            "liftOver_needed": "n/a",
            "notes": "Use only when mapping to hg19 resources",
        },
        {
            "dataset": "BrainSeq WGBS meQTL (external)",
            "path_template": "external; see public_meqtl_resources.tsv",
            "genome_build": "hg38 (confirm on download)",
            "coordinate_system": "CpG / SNP positions",
            "evidence": "Perzel Mandell et al. 2021 Nat Commun",
            "liftOver_needed": "maybe",
            "notes": "Primary external resource; verify build on ingest",
        },
        {
            "dataset": "Jaffe DLPFC 450K meQTL (external)",
            "path_template": "external; see public_meqtl_resources.tsv",
            "genome_build": "hg19 likely",
            "coordinate_system": "Illumina cg IDs + positions",
            "evidence": "published prefrontal meQTL catalog",
            "liftOver_needed": "yes_if_hg19",
            "notes": "Harmonize cg IDs to WGBS CpGs via genomic position",
        },
        {
            "dataset": "Hippocampus array meQTL (Schulz 2017)",
            "path_template": "external; see public_meqtl_resources.tsv",
            "genome_build": "hg19 likely",
            "coordinate_system": "array CpG positions",
            "evidence": "Nat Commun 2017",
            "liftOver_needed": "yes_if_hg19",
            "notes": "",
        },
        {
            "dataset": "S-LDSC schizophrenia sumstats (existing)",
            "path_template": "heritability/.../clinical_enrichment/s-ldsc/hg19/",
            "genome_build": "hg19",
            "coordinate_system": "GWAS variants",
            "evidence": "existing pipeline under hg19/",
            "liftOver_needed": "already handled in existing pipeline",
            "notes": "Not primary AJHG strengthening path",
        },
    ]
    # existence flags
    for row in rows:
        tmpl = row["path_template"]
        if tmpl.startswith("external") or tmpl.startswith("heritability/..."):
            row["example_exists"] = "n/a"
            continue
        exists_any = False
        if "{region}" in tmpl:
            for region in cfg["regions"]:
                path = project / tmpl.format(region=region, chrom=1)
                if not path.exists() and "cpg" in tmpl:
                    path = Path(cfg["cpg_methylation_root"]) / tmpl.format(region=region, chrom=1)
                if path.exists():
                    exists_any = True
                    break
        else:
            path = project / tmpl
            exists_any = path.exists()
        row["example_exists"] = str(exists_any)

    write_tsv(
        OUTDIR / "genome_build_audit.tsv",
        rows,
        [
            "dataset", "path_template", "genome_build", "coordinate_system",
            "evidence", "liftOver_needed", "example_exists", "notes",
        ],
    )


def write_public_meqtl_resources() -> None:
    rows = [
        {
            "resource_id": "brainseq_wgbs_meqtl",
            "citation": "Perzel Mandell et al. Nat Commun 2021; PMID:34471112",
            "tissue_regions": "DLPFC,hippocampus",
            "methylation_platform": "WGBS",
            "sample_size": "DLPFC n=165; hippocampus n=179; 183 donors",
            "ancestry_notes": "LIBD adult postmortem; ancestry composition documented in paper",
            "genome_build": "verify_on_download",
            "cis_window": "as published",
            "significance_threshold": "as published (cis-meQTL catalogs)",
            "access_url": "https://eqtl.brainseq.org/WGBS_meQTL/",
            "data_download": "PsychENCODE / BrainSeq portals; browser at eqtl.brainseq.org",
            "overlap_with_vmr_cpgs": "expected_high_platform_matched",
            "priority": "primary",
            "inclusion_decision": "include",
            "harmonization_notes": "Prefer genomic coordinates over probe IDs; strand/allele for effect-direction comparisons",
        },
        {
            "resource_id": "jaffe_dlpfc_450k_meqtl",
            "citation": "Jaffe et al. Nat Neurosci 2016 / related LIBD DLPFC 450K meQTL releases",
            "tissue_regions": "DLPFC",
            "methylation_platform": "Illumina_450K",
            "sample_size": "adult DLPFC ~258 non-psychiatric in cited catalogs",
            "ancestry_notes": "predominantly European-ancestry LIBD cohorts in published catalogs",
            "genome_build": "hg19_likely",
            "cis_window": "as published",
            "significance_threshold": "FDR 1% in cited catalogs",
            "access_url": "https://www.nature.com/articles/nn.4181",
            "data_download": "supplement / LIBD resources; map cg IDs to hg38 positions",
            "overlap_with_vmr_cpgs": "moderate_array_subset",
            "priority": "secondary",
            "inclusion_decision": "include",
            "harmonization_notes": "liftOver positions if hg19; do not pool with WGBS without justification",
        },
        {
            "resource_id": "schulz_hippocampus_array_meqtl",
            "citation": "Schulz et al. Nat Commun 2017; DOI 10.1038/s41467-017-01818-4",
            "tissue_regions": "hippocampus",
            "methylation_platform": "Illumina_array",
            "sample_size": "as published",
            "ancestry_notes": "as published",
            "genome_build": "hg19_likely",
            "cis_window": "as published",
            "significance_threshold": "FDR as published (~14k cis-meQTL CpGs)",
            "access_url": "https://doi.org/10.1038/s41467-017-01818-4",
            "data_download": "Supplementary Data / sciebo link in paper",
            "overlap_with_vmr_cpgs": "moderate_array_subset",
            "priority": "secondary",
            "inclusion_decision": "include",
            "harmonization_notes": "liftOver if needed; hippocampus-specific external support",
        },
        {
            "resource_id": "hannon_fetal_brain_meqtl",
            "citation": "Hannon et al. fetal brain meQTL (cited in Schulz 2017)",
            "tissue_regions": "fetal_brain",
            "methylation_platform": "Illumina_450K",
            "sample_size": "166 fetal brains (cited)",
            "ancestry_notes": "as published",
            "genome_build": "hg19_likely",
            "cis_window": "as published",
            "significance_threshold": "Bonferroni as published",
            "access_url": "",
            "data_download": "as published supplements",
            "overlap_with_vmr_cpgs": "limited_developmental_context",
            "priority": "exploratory",
            "inclusion_decision": "optional_positive_control_only",
            "harmonization_notes": "Not adult postmortem; do not use as primary external validation",
        },
    ]
    write_tsv(
        OUTDIR / "public_meqtl_resources.tsv",
        rows,
        [
            "resource_id", "citation", "tissue_regions", "methylation_platform",
            "sample_size", "ancestry_notes", "genome_build", "cis_window",
            "significance_threshold", "access_url", "data_download",
            "overlap_with_vmr_cpgs", "priority", "inclusion_decision",
            "harmonization_notes",
        ],
    )


def write_inventory_summaries(cfg: dict) -> None:
    project = Path(cfg["project_root"])
    staff = Path(cfg["cpg_methylation_root"])
    rows = []
    for region in cfg["regions"]:
        vmr = project / cfg["vmr_bed_template"].format(region=region)
        pred = project / cfg["local_predictability_summary_template"].format(region=region)
        cpg_dir = staff / "vmr-analysis" / region / "_m" / "cpg"
        n_chr = len(list(cpg_dir.glob("chr_*"))) if cpg_dir.exists() else 0
        rows.append({
            "region": region,
            "n_vmrs": count_lines(vmr) if vmr.exists() else "",
            "n_predictability_rows": (count_lines(pred) - 1) if pred.exists() and isinstance(count_lines(pred), int) else "",
            "n_staff_cpg_chr_dirs": n_chr,
            "predictability_path": str(pred),
            "vmr_path": str(vmr),
            "staff_cpg_dir": str(cpg_dir),
        })
    write_tsv(
        OUTDIR / "region_input_summary.tsv",
        rows,
        ["region", "n_vmrs", "n_predictability_rows", "n_staff_cpg_chr_dirs", "predictability_path", "vmr_path", "staff_cpg_dir"],
    )


def write_reproducibility() -> None:
    git_commit = ""
    try:
        import subprocess
        git_commit = subprocess.check_output(
            ["git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        git_commit = "unavailable"
    write_tsv(
        OUTDIR / "data_dictionary_reproducibility.tsv",
        [{
            "script": str(Path(__file__).resolve()),
            "execution_utc": datetime.now(timezone.utc).isoformat(),
            "python_version": sys.version.replace("\n", " "),
            "platform": platform.platform(),
            "git_commit": git_commit,
            "paths_config": str(PATHS_YML),
            "pyyaml_available": str(yaml is not None),
        }],
        ["script", "execution_utc", "python_version", "platform", "git_commit", "paths_config", "pyyaml_available"],
    )


def main() -> None:
    cfg = load_paths()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_sample_overlap(cfg)
    write_covariate_dictionary(cfg)
    write_genome_build_audit(cfg)
    write_public_meqtl_resources()
    write_inventory_summaries(cfg)
    write_reproducibility()
    print(f"Wrote data dictionary tables to {OUTDIR}")


if __name__ == "__main__":
    main()
