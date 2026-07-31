#!/usr/bin/env python3
"""Shared helpers for meQTL-validation analysis modules."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = PROJECT_ROOT / "config"


def load_yaml(name: str) -> dict[str, Any]:
    path = CONFIG_DIR / name
    if yaml is None:
        raise SystemExit(f"PyYAML required to load {path}")
    with path.open() as handle:
        return yaml.safe_load(handle)


def load_paths() -> dict[str, Any]:
    return load_yaml("paths.yml")


def cohort_key(population: str) -> str:
    """Map analysis population to path cohort: AA vs all_individuals (EA/expanded)."""
    pop = (population or "AA").upper()
    if pop == "AA":
        return "AA"
    if pop in {"EA", "ALL", "ALL_INDIVIDUALS", "EXPANDED"}:
        return "all_individuals"
    raise ValueError(f"Unknown population/cohort: {population}")


def genotype_paths(paths: dict[str, Any], population: str = "AA") -> dict[str, str]:
    """Return genotype file dict for AA or all_individuals/EA."""
    key = "AA" if cohort_key(population) == "AA" else "EA"
    geno = paths.get("genotype", {})
    if key in geno:
        return geno[key]
    if key == "EA" and "all_individuals" in geno:
        return geno["all_individuals"]
    raise KeyError(f"genotype.{key} missing from paths.yml")


def cpg_matrix_relpath(paths: dict[str, Any], region: str, chrom: str, population: str = "AA") -> str:
    """Relative path (under cpg_methylation_root) to non-residualized cpg_meth.phen."""
    chrom = str(chrom).replace("chr", "")
    if cohort_key(population) == "AA":
        tmpl = paths.get(
            "cpg_methylation_matrix_aa_template",
            paths.get("cpg_methylation_matrix_template"),
        )
    else:
        tmpl = paths.get(
            "cpg_methylation_matrix_ea_template",
            paths.get("cpg_methylation_matrix_all_individuals_template"),
        )
    if not tmpl:
        raise KeyError("CpG methylation matrix template missing from paths.yml")
    return tmpl.format(region=region, chrom=chrom)


def cpg_samples_relpath(paths: dict[str, Any], region: str, population: str = "AA") -> str:
    if cohort_key(population) == "AA":
        tmpl = paths.get("cpg_samples_aa_template", f"vmr-analysis/{region}/_m/samples.txt")
    else:
        tmpl = paths.get(
            "cpg_samples_all_individuals_template",
            f"vmr-analysis/all_individuals/{region}/_m/samples.txt",
        )
    return tmpl.format(region=region)


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def norm_region(value: str) -> str:
    lower = (value or "").strip().lower()
    if lower == "caudate nucleus":
        return "caudate"
    return lower


def nonmissing(row: dict, col: str) -> bool:
    return row.get(col, "") not in ("", "NA", "NaN", "nan", "NULL", "None")


def read_psam_brnums(psam: Path) -> set[str]:
    ids: set[str] = set()
    with psam.open() as handle:
        first = handle.readline().rstrip("\n").split("\t")
        headerish = [h.lstrip("#") for h in first]
        has_header = "FID" in headerish or "IID" in headerish
        if has_header:
            fid_idx = headerish.index("FID") if "FID" in headerish else 0
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) > fid_idx:
                    ids.add(parts[fid_idx])
        else:
            ids.add(first[0])
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if parts:
                    ids.add(parts[0])
    return ids


def chrom_label(chrom: str | int) -> str:
    c = str(chrom).replace("chr", "")
    return f"chr{c}"
