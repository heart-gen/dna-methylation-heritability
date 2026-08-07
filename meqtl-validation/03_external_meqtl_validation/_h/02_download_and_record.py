#!/usr/bin/env python3
"""Record downloaded external meQTL archives into download_manifest.tsv.

Does not overwrite existing checksum rows for the same filename.
BrainSeq full cis catalogs require Synapse (syn25992404); this script only
records files already present under _m/raw/.
"""

from __future__ import annotations

import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import PROJECT_ROOT, read_tsv, write_tsv  # noqa: E402

OUTDIR = PROJECT_ROOT / "meqtl-validation" / "03_external_meqtl_validation" / "_m"
RAW = OUTDIR / "raw"

# Known provenance for files already staged under raw/
KNOWN_URLS = {
    "GSE74193_jaffe_onlineTable2_meQTLs_allPairs.csv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74193/suppl/"
        "GSE74193_jaffe_onlineTable2_meQTLs_allPairs.csv.gz"
    ),
    "hg19ToHg38.over.chain.gz": (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
    ),
}


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            buf = handle.read(chunk)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()


def infer_resource_id(path: Path) -> str:
    try:
        rel = path.relative_to(RAW)
        return rel.parts[0]
    except ValueError:
        return "support"


def infer_build(path: Path) -> str:
    name = path.name.lower()
    if "hg19" in name or "450k" in name or "jaffe" in name or "schulz" in name or "gse74193" in name:
        return "hg19"
    if "brainseq" in str(path).lower() or "41467_2021" in name:
        return "hg38_likely"
    return "unknown"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    manifest_path = OUTDIR / "download_manifest.tsv"
    existing = read_tsv(manifest_path) if manifest_path.exists() else []
    seen = {(r.get("resource_id"), r.get("filename")) for r in existing}

    rows = list(existing)
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    candidates: list[Path] = []
    for rid_dir in sorted(RAW.glob("*")):
        if not rid_dir.is_dir():
            continue
        for path in rid_dir.rglob("*"):
            if path.is_file() and path.stat().st_size > 0:
                candidates.append(path)
    support = OUTDIR / "support"
    if support.exists():
        candidates.extend(p for p in support.glob("*") if p.is_file())

    for path in candidates:
        rid = infer_resource_id(path) if path.is_relative_to(RAW) else "support"
        key = (rid, path.name)
        if key in seen:
            continue
        url = KNOWN_URLS.get(path.name, "")
        if not url and "MOESM" in path.name:
            url = "nature_supplement"
        rows.append({
            "resource_id": rid,
            "filename": path.name,
            "url": url or str(path),
            "download_utc": now,
            "sha256": sha256_file(path),
            "genome_build_as_downloaded": infer_build(path),
            "notes": f"bytes={path.stat().st_size}; path={path}",
        })
        seen.add(key)

    write_tsv(
        manifest_path,
        rows,
        ["resource_id", "filename", "url", "download_utc", "sha256", "genome_build_as_downloaded", "notes"],
    )
    print(f"Wrote {len(rows)} manifest rows to {manifest_path}")


if __name__ == "__main__":
    main()
