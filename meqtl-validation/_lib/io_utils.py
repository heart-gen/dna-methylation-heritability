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
