#!/usr/bin/env python3
"""Audit available inputs for the follow-up pilot.

This script intentionally does not fit phenotype-wide models. It records
available files, sample overlap, missingness, and high-level feature counts so
later analyses can fail early if an expected input is absent or inconsistent.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import platform
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[3]
CONFIG = ANALYSIS_DIR / "_h" / "analysis_inputs.json"
OUTDIR = ANALYSIS_DIR / "_m"


def load_config() -> dict:
    with CONFIG.open() as handle:
        return json.load(handle)


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


def exists_info(root_name: str, root: Path, rel_path: str) -> dict:
    path = root / rel_path
    info = {
        "root": root_name,
        "relative_path": rel_path,
        "path": str(path),
        "exists": os.path.exists(path),
        "is_file": path.is_file(),
        "is_dir": path.is_dir(),
        "is_symlink": path.is_symlink(),
        "symlink_target": os.readlink(path) if path.is_symlink() else "",
        "size_bytes": "",
    }
    if info["exists"] and path.is_file():
        try:
            info["size_bytes"] = path.stat().st_size
        except OSError:
            info["size_bytes"] = ""
    return info


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


def norm_region(value: str) -> str:
    lower = (value or "").strip().lower()
    if lower == "caudate nucleus":
        return "caudate"
    return lower


def nonmissing(row: dict, col: str) -> bool:
    val = row.get(col, "")
    return val not in ("", "NA", "NaN", "nan", "NULL", "None")


def audit_inventory(cfg: dict) -> None:
    project = Path(cfg["project_root"])
    staff = Path(cfg["staff_project_root"])
    roots = [("project", project), ("staff", staff)]

    rels: list[str] = [
        cfg["phenotype_table"],
        "inputs/phenotypes/_m/phenotypes-AA.tsv",
        "inputs/phenotypes/_m/phenotypes-DNAm.tsv",
        "inputs/phenotypes/_m/phenotypes-DNAm-all.tsv",
        "inputs/phenotypes/_m/phenotypes-all.tsv",
        "inputs/supportfiles/_m/repeat-masker-hg38.gz",
        "inputs/supportfiles/_m/hg38ToHg19.over.chain",
    ]
    rels.extend(cfg["region_phenotype_tables"].values())
    rels.extend(cfg["genotype_files"].values())

    for region in cfg["regions"]:
        rels.append(cfg["vmr_bed_template"].format(region=region))
        rels.append(cfg["local_predictability_summary_template"].format(region=region))
        rels.append(cfg["cell_proportions_template"].format(region=region))
        rels.append(f"inputs/wgbs-data/combined/_m/{region}_bsseq_h5se")
        rels.append(f"vmr-analysis/{region}/_m/{region}_cpg_meth.tar.gz")
        rels.append(f"vmr-analysis/{region}/_m/cpg/top1_cpg.tsv")

    rows = []
    for root_name, root in roots:
        for rel in sorted(set(rels)):
            rows.append(exists_info(root_name, root, rel))

    write_tsv(
        OUTDIR / "input_inventory.tsv",
        rows,
        ["root", "relative_path", "path", "exists", "is_file", "is_dir", "is_symlink", "symlink_target", "size_bytes"],
    )


def audit_phenotypes(cfg: dict) -> None:
    phen_path = Path(cfg["project_root"]) / cfg["phenotype_table"]
    rows = read_table(phen_path)
    for row in rows:
        row["_region_norm"] = norm_region(row.get("region", ""))

    count_rows = []
    counter = Counter((row.get("race", ""), row["_region_norm"], row.get("primarydx", "")) for row in rows)
    for (race, region, dx), n in sorted(counter.items()):
        count_rows.append({"race": race, "region": region, "primarydx": dx, "n_samples": n})
    write_tsv(OUTDIR / "phenotype_counts.tsv", count_rows, ["race", "region", "primarydx", "n_samples"])

    donor_region = defaultdict(set)
    for row in rows:
        donor_region[row.get("brnum", "")].add(row["_region_norm"])
    repeated_rows = []
    for n_regions, n_donors in sorted(Counter(len(v) for v in donor_region.values()).items()):
        repeated_rows.append({"n_regions_observed": n_regions, "n_donors": n_donors})
    write_tsv(OUTDIR / "repeated_donor_counts.tsv", repeated_rows, ["n_regions_observed", "n_donors"])

    covariates = [
        "agedeath", "sex", "primarydx", "pmi", "ph",
        "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5",
        "snpPC6", "snpPC7", "snpPC8", "snpPC9", "snpPC10",
        "smoking", "nicotine", "cocaine", "ethanol", "amphetamines",
        "codeine", "morphine", "fentanyl", "antipsychotics",
        "lifetime_antipsych", "education", "manner_of_death",
    ]
    miss_rows = []
    for race in sorted({row.get("race", "") for row in rows}):
        for region in sorted({row["_region_norm"] for row in rows}):
            subset = [row for row in rows if row.get("race", "") == race and row["_region_norm"] == region]
            if not subset:
                continue
            for covar in covariates:
                if covar not in rows[0]:
                    continue
                missing = sum(1 for row in subset if not nonmissing(row, covar))
                miss_rows.append({
                    "race": race,
                    "region": region,
                    "covariate": covar,
                    "n_samples": len(subset),
                    "n_missing": missing,
                    "pct_missing": round(100 * missing / len(subset), 3),
                })
    write_tsv(OUTDIR / "covariate_missingness.tsv", miss_rows, ["race", "region", "covariate", "n_samples", "n_missing", "pct_missing"])

    primary_covars = ["agedeath", "sex", "primarydx", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5"]
    complete_rows = []
    for race in sorted({row.get("race", "") for row in rows}):
        for region in sorted({row["_region_norm"] for row in rows}):
            subset = [row for row in rows if row.get("race", "") == race and row["_region_norm"] == region]
            if not subset:
                continue
            complete = [row for row in subset if all(nonmissing(row, covar) for covar in primary_covars)]
            dx_counts = Counter(row.get("primarydx", "") for row in complete)
            complete_rows.append({
                "race": race,
                "region": region,
                "analysis": "primary_schizophrenia_or_aging_covariate_set",
                "required_covariates": ",".join(primary_covars),
                "n_total": len(subset),
                "n_complete": len(complete),
                "n_control_complete": dx_counts.get("Control", 0),
                "n_schizophrenia_complete": dx_counts.get("Schizo", 0),
            })
    write_tsv(
        OUTDIR / "primary_model_complete_cases.tsv",
        complete_rows,
        ["race", "region", "analysis", "required_covariates", "n_total", "n_complete", "n_control_complete", "n_schizophrenia_complete"],
    )


def audit_region_files(cfg: dict) -> None:
    project = Path(cfg["project_root"])
    staff = Path(cfg["staff_project_root"])
    roots = [("project", project), ("staff", staff)]

    vmr_rows = []
    pred_rows = []
    meth_rows = []
    cell_rows = []

    for root_name, root in roots:
        for region in cfg["regions"]:
            vmr_rel = cfg["vmr_bed_template"].format(region=region)
            vmr_path = root / vmr_rel
            vmr_rows.append({
                "root": root_name,
                "region": region,
                "path": str(vmr_path),
                "exists": vmr_path.exists(),
                "n_vmrs": count_lines(vmr_path) if vmr_path.exists() else "",
            })

            pred_rel = cfg["local_predictability_summary_template"].format(region=region)
            pred_path = root / pred_rel
            pred_lines = count_lines(pred_path) if pred_path.exists() else ""
            pred_rows.append({
                "root": root_name,
                "region": region,
                "path": str(pred_path),
                "exists": pred_path.exists(),
                "n_rows_excluding_header": pred_lines - 1 if isinstance(pred_lines, int) else "",
                "n_columns": header_n_columns(pred_path) if pred_path.exists() else "",
            })

            cell_rel = cfg["cell_proportions_template"].format(region=region)
            cell_path = root / cell_rel
            if cell_path.exists():
                rows = read_table(cell_path)
                cell_rows.append({
                    "root": root_name,
                    "region": region,
                    "path": str(cell_path),
                    "exists": True,
                    "n_rows": len(rows),
                    "n_samples": len({row.get("sample_id", "") for row in rows}),
                    "cell_types": ",".join(sorted({row.get("cell_type", "") for row in rows})),
                })
            else:
                cell_rows.append({"root": root_name, "region": region, "path": str(cell_path), "exists": False, "n_rows": "", "n_samples": "", "cell_types": ""})

            cpg_dir = root / "vmr-analysis" / region / "_m" / "cpg"
            chr_dirs = sorted(cpg_dir.glob("chr_*")) if cpg_dir.exists() else []
            cpg_files = [path / "cpg_meth.phen" for path in chr_dirs if (path / "cpg_meth.phen").exists()]
            res_files = [path / "res_cpg_meth.phen" for path in chr_dirs if (path / "res_cpg_meth.phen").exists()]
            example = cpg_files[0] if cpg_files else None
            meth_rows.append({
                "root": root_name,
                "region": region,
                "cpg_dir": str(cpg_dir),
                "cpg_dir_exists": cpg_dir.exists(),
                "n_chr_dirs": len(chr_dirs),
                "n_cpg_meth_phen": len(cpg_files),
                "n_residualized_cpg_meth_phen": len(res_files),
                "example_cpg_meth_path": str(example) if example else "",
                "example_cpg_columns_including_FID_IID": header_n_columns(example) if example else "",
                "archive_exists": (root / "vmr-analysis" / region / "_m" / f"{region}_cpg_meth.tar.gz").exists(),
            })

    write_tsv(OUTDIR / "vmr_inventory.tsv", vmr_rows, ["root", "region", "path", "exists", "n_vmrs"])
    write_tsv(OUTDIR / "local_predictability_inventory.tsv", pred_rows, ["root", "region", "path", "exists", "n_rows_excluding_header", "n_columns"])
    write_tsv(OUTDIR / "cell_composition_inventory.tsv", cell_rows, ["root", "region", "path", "exists", "n_rows", "n_samples", "cell_types"])
    write_tsv(
        OUTDIR / "cpg_methylation_inventory.tsv",
        meth_rows,
        ["root", "region", "cpg_dir", "cpg_dir_exists", "n_chr_dirs", "n_cpg_meth_phen", "n_residualized_cpg_meth_phen", "example_cpg_meth_path", "example_cpg_columns_including_FID_IID", "archive_exists"],
    )


def audit_genotypes(cfg: dict) -> None:
    project = Path(cfg["project_root"])
    rows = []
    for label, rel in cfg["genotype_files"].items():
        path = project / rel
        row = exists_info("project", project, rel)
        row["label"] = label
        row["n_rows_excluding_header"] = ""
        row["n_columns"] = ""
        if path.exists() and path.is_file() and path.suffix in {".psam", ".eigenvec", ".eigenval"}:
            lines = count_lines(path)
            row["n_rows_excluding_header"] = lines - 1 if isinstance(lines, int) and path.suffix == ".psam" else lines
            row["n_columns"] = header_n_columns(path)
        rows.append(row)
    write_tsv(
        OUTDIR / "genotype_inventory.tsv",
        rows,
        ["label", "root", "relative_path", "path", "exists", "is_file", "is_dir", "is_symlink", "symlink_target", "size_bytes", "n_rows_excluding_header", "n_columns"],
    )


def audit_molecular_links(cfg: dict) -> None:
    base = Path(cfg["project_root"]) / cfg["regulatory_context_root"]
    rows = []
    pattern = str(base / "**" / "sample_feature_coverage.tsv")
    for path_str in sorted(glob.glob(pattern, recursive=True)):
        path = Path(path_str)
        rel_parts = path.relative_to(base).parts
        rows_in = read_table(path)
        first = rows_in[0] if rows_in else {}
        rows.append({
            "path": str(path),
            "relative_context": str(path.relative_to(base)),
            "n_rows": len(rows_in),
            "cohort": first.get("cohort", ""),
            "tissue": first.get("tissue", rel_parts[0] if rel_parts else ""),
            "modality": first.get("modality", ""),
            "population": first.get("population", ""),
            "n_control_pheno": first.get("n_control_pheno", ""),
            "n_rna_samples": first.get("n_rna_samples", ""),
            "n_model_samples_pre_pair_filter": first.get("n_model_samples_pre_pair_filter", ""),
            "n_features_linked": first.get("n_features_linked", ""),
            "model": first.get("model", ""),
        })
    write_tsv(
        OUTDIR / "molecular_link_coverage.tsv",
        rows,
        ["path", "relative_context", "n_rows", "cohort", "tissue", "modality", "population", "n_control_pheno", "n_rna_samples", "n_model_samples_pre_pair_filter", "n_features_linked", "model"],
    )


def audit_reproducibility() -> None:
    rows = [{
        "script": str(Path(__file__).resolve()),
        "execution_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "random_seed": "not used; deterministic file audit",
        "git_commit": "not recorded; git status failed under sandbox because Git LFS requires writing .git/lfs/tmp",
        "config": str(CONFIG),
    }]
    write_tsv(OUTDIR / "audit_reproducibility.tsv", rows, ["script", "execution_utc", "python_version", "platform", "random_seed", "git_commit", "config"])


def main() -> None:
    cfg = load_config()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    audit_inventory(cfg)
    audit_phenotypes(cfg)
    audit_region_files(cfg)
    audit_genotypes(cfg)
    audit_molecular_links(cfg)
    audit_reproducibility()


if __name__ == "__main__":
    main()
