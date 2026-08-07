#!/usr/bin/env python3
"""Initialize external meQTL validation workspace and overlap scaffold.

Does not download remote resources automatically (manual/authenticated portals).
Writes a download checklist from the Phase 0 catalog.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import PROJECT_ROOT, read_tsv, write_tsv  # noqa: E402


def main() -> None:
    catalog = PROJECT_ROOT / "inputs" / "data_dictionary" / "_m" / "public_meqtl_resources.tsv"
    outdir = PROJECT_ROOT / "meqtl-validation" / "03_external_meqtl_validation" / "_m"
    raw = outdir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(catalog)
    checklist = []
    for r in rows:
        rid = r["resource_id"]
        dest = raw / rid
        dest.mkdir(parents=True, exist_ok=True)
        checklist.append({
            "resource_id": rid,
            "priority": r.get("priority", ""),
            "inclusion_decision": r.get("inclusion_decision", ""),
            "access_url": r.get("access_url", ""),
            "local_raw_dir": str(dest),
            "download_complete": "false",
            "harmonized_path": "",
            "notes": "Place original files here; record checksum/date in download_manifest.tsv",
        })
    write_tsv(outdir / "download_checklist.tsv", checklist)
    if not (outdir / "download_manifest.tsv").exists():
        write_tsv(
            outdir / "download_manifest.tsv",
            [],
            ["resource_id", "filename", "url", "download_utc", "sha256", "genome_build_as_downloaded", "notes"],
        )
    print(f"Initialized external meQTL workspace under {outdir}")


if __name__ == "__main__":
    main()
