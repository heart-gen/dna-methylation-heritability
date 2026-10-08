#!/usr/bin/env python3
"""Stage the Rizzardi et al. 2019 neuronal CG-DMRs as hg19 BED assets.

Usage:
    python -I inputs/supportfiles/_h/03_build_rizzardi_dmr_asset.py --download-dir <empty dir>

Reads `config/sldsc_external_annotations.yml`. For every annotation with a
`source_url` ending in .bb, downloads it into --download-dir, checks the SHA-256
the config pins, converts the bigBed to a three-column BED with pyBigWig, checks
the interval count against `expected_n_intervals`, and writes the `bed` path the
config names. Annotations without a downloadable source (status pending_access)
are skipped and reported.

The source is the authors' own UCSC hub (hg19), which is the build the LDSC
reference panels use, so no liftover is involved. Intervals are sorted and
merged so that membership is a single overlap test downstream.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import subprocess
import sys
import urllib.request
from pathlib import Path

import pyBigWig
import yaml


def repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").is_dir():
            return parent
    raise SystemExit("Could not locate repository root")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--download-dir", required=True, type=Path)
    args = ap.parse_args()
    root = repo_root()
    cfg = yaml.safe_load((root / "config" / "sldsc_external_annotations.yml").read_text())
    args.download_dir.mkdir(parents=True, exist_ok=True)

    report = [("annotation", "status", "n_intervals_source", "n_autosomal",
               "n_merged", "bp_merged", "bed_sha256")]
    for a in cfg["annotations"]:
        url = a.get("source_url") or ""
        if not url.endswith(".bb") or not a.get("bed"):
            report.append((a["name"], a.get("status", "no_downloadable_source"),
                           "", "", "", "", ""))
            continue
        src = args.download_dir / Path(url).name
        if not src.exists():
            urllib.request.urlretrieve(url, src)
        got = sha256(src)
        if got != a["source_sha256"]:
            raise SystemExit(f"{a['name']}: {src} sha256 {got} != pinned {a['source_sha256']}")
        bb = pyBigWig.open(str(src))
        if not bb.isBigBed():
            raise SystemExit(f"{src} is not a bigBed")
        rows = []
        for chrom, length in bb.chroms().items():
            for s, e, _ in bb.entries(chrom, 0, length) or []:
                rows.append((chrom, s, e))
        if len(rows) != int(a["expected_n_intervals"]):
            raise SystemExit(f"{a['name']}: {len(rows)} intervals, config expects "
                             f"{a['expected_n_intervals']}")
        auto = sorted((r for r in rows if r[0] in AUTOSOMES),
                      key=lambda r: (int(r[0][3:]), r[1], r[2]))
        merged = []
        for c, s, e in auto:
            if merged and merged[-1][0] == c and s <= merged[-1][2]:
                merged[-1][2] = max(merged[-1][2], e)
            else:
                merged.append([c, s, e])
        out = root / a["bed"]
        # mtime=0 keeps the gzip header fixed, so a rebuild is byte-identical
        # and its SHA-256 can be pinned.
        with open(out, "wb") as raw, \
                gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
            for c, s, e in merged:
                gz.write(f"{c}\t{s}\t{e}\n".encode())
        report.append((a["name"], "ready", len(rows), len(auto), len(merged),
                       sum(e - s for _, s, e in merged), sha256(out)))

    rep = root / "inputs" / "supportfiles" / "_m" / "rizzardi2019-cgdmr-asset-report.tsv"
    with open(rep, "w") as fh:
        for r in report:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(open(rep).read(), file=sys.stderr)


if __name__ == "__main__":
    main()
