#!/usr/bin/env python3
"""Stage Phase 6 technical-annotation assets under inputs/supportfiles/_m.

Downloads / converts (hg38):
  - ENCODE blacklist v2
  - UCSC genomicSuperDups
  - Hoffman k24 Umap multi-track mappability bigWig
  - LIBD C/T SNP positions at CpGs (DEM2) → BED
  - RepeatMasker LINE/L1 subset (from existing repeat-masker-hg38.gz)
  - Common autosomal SNP ±snp_proximity_bp windows from AA pvar (MAF filter)

Writes asset_manifest.tsv with checksums and notes. Does not overwrite
existing non-empty files unless --force.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import PROJECT_ROOT, load_paths, load_yaml, write_tsv  # noqa: E402

SUPPORT = PROJECT_ROOT / "inputs" / "supportfiles" / "_m"
CT_SNP_DIR = Path("/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--force", action="store_true")
    p.add_argument("--skip-umap", action="store_true", help="Skip large mappability bigWig download")
    p.add_argument("--skip-snp-windows", action="store_true")
    return p.parse_args()


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            buf = handle.read(chunk)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()


def download(url: str, dest: Path, force: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"exists: {dest}")
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"download: {url} -> {dest}")
    subprocess.check_call(["curl", "-L", "--fail", "-o", str(tmp), url])
    tmp.replace(dest)


def convert_segdups(src: Path, dest: Path, force: bool) -> None:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"exists: {dest}")
        return
    print(f"convert segdups -> {dest}")
    # UCSC genomicSuperDups.txt.gz columns:
    # bin chrom chromStart chromEnd name ...
    n = 0
    with gzip.open(src, "rt") as hin, gzip.open(dest, "wt") as hout:
        for line in hin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end = parts[1], parts[2], parts[3]
            if not str(chrom).startswith("chr"):
                continue
            hout.write(f"{chrom}\t{start}\t{end}\tsegdup\n")
            n += 1
    print(f"  wrote {n:,} segdup intervals")


def convert_ct_snps(dest: Path, force: bool) -> None:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"exists: {dest}")
        return
    print(f"convert CT SNPs -> {dest}")
    with gzip.open(dest, "wt") as hout:
        for chrom in range(1, 23):
            src = CT_SNP_DIR / f"chr{chrom}"
            if not src.exists():
                print(f"WARNING: missing {src}")
                continue
            with src.open() as hin:
                for line in hin:
                    pos = line.strip()
                    if not pos:
                        continue
                    p = int(float(pos))
                    # 0-based half-open for CpG C position
                    hout.write(f"chr{chrom}\t{p-1}\t{p}\tCT_SNP\n")


def extract_line_l1(repeat_gz: Path, dest: Path, force: bool) -> None:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"exists: {dest}")
        return
    print(f"extract LINE/L1 from {repeat_gz}")
    with gzip.open(repeat_gz, "rt") as hin, gzip.open(dest, "wt") as hout:
        for line in hin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            name = parts[3]
            if name.startswith("L1") or name.startswith("LINE"):
                hout.write(line if line.endswith("\n") else line + "\n")


def build_snp_proximity_windows(
    pvar: Path,
    dest: Path,
    maf_min: float,
    window_bp: int,
    force: bool,
) -> None:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"exists: {dest}")
        return
    print(f"build SNP ±{window_bp}bp windows from {pvar} (MAF>={maf_min})")
    # pvar may lack AF; use plink2 if available to get AF, else keep all autosomal SNPs
    # Prefer existing filtered meqtl genotypes AF if present
    n_written = 0
    open_pvar = pvar.open
    with open_pvar() as hin, gzip.open(dest, "wt") as hout:
        header = None
        af_idx = None
        for line in hin:
            if line.startswith("##"):
                continue
            if line.startswith("#") or header is None:
                header = line.lstrip("#").rstrip("\n").split("\t")
                for cand in ("AF", "ALT_FREQS", "MAF"):
                    if cand in header:
                        af_idx = header.index(cand)
                        break
                # INFO may contain AF=
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            if chrom.startswith("chr"):
                cnorm = chrom
            else:
                if chrom not in {str(i) for i in range(1, 23)}:
                    continue
                cnorm = f"chr{chrom}"
            if cnorm.replace("chr", "") not in {str(i) for i in range(1, 23)}:
                continue
            try:
                pos = int(parts[1])
            except (ValueError, IndexError):
                continue
            if af_idx is not None:
                try:
                    af = float(parts[af_idx].split(",")[0])
                    maf = min(af, 1.0 - af)
                    if maf < maf_min:
                        continue
                except ValueError:
                    pass
            elif len(parts) > 7 and "AF=" in parts[7]:
                try:
                    info = parts[7]
                    af = float(info.split("AF=")[1].split(";")[0].split(",")[0])
                    maf = min(af, 1.0 - af)
                    if maf < maf_min:
                        continue
                except (IndexError, ValueError):
                    pass
            start = max(0, pos - 1 - window_bp)
            end = pos + window_bp
            hout.write(f"{cnorm}\t{start}\t{end}\tcommon_snp_window\n")
            n_written += 1
            if n_written % 2_000_000 == 0:
                print(f"  ... {n_written:,} windows")
    print(f"wrote {n_written:,} SNP proximity windows")


def main() -> None:
    args = parse_args()
    SUPPORT.mkdir(parents=True, exist_ok=True)
    paths = load_paths()
    meqtl = load_yaml("meqtl_parameters.yml")
    thresholds = load_yaml("analysis_thresholds.yml")
    window_bp = int(thresholds["repeat_sensitivity"]["snp_proximity_bp"])
    maf_min = float(meqtl["genotype_qc"]["maf_min"])
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    blacklist = SUPPORT / "hg38-blacklist.v2.bed.gz"
    download(
        "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
        blacklist,
        args.force,
    )

    segdup_txt = SUPPORT / "genomicSuperDups.txt.gz"
    download(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz",
        segdup_txt,
        args.force,
    )
    segdup_bed = SUPPORT / "genomicSuperDups.hg38.bed.gz"
    convert_segdups(segdup_txt, segdup_bed, args.force)

    umap_bw = SUPPORT / "k24.Umap.MultiTrackMappability.bw"
    if not args.skip_umap:
        download(
            "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k24.Umap.MultiTrackMappability.bw",
            umap_bw,
            args.force,
        )

    ct_bed = SUPPORT / "libd_CT_snps_at_cpg.hg38.bed.gz"
    convert_ct_snps(ct_bed, args.force)

    repeat_gz = SUPPORT / "repeat-masker-hg38.gz"
    line_l1 = SUPPORT / "repeat-masker.LINE_L1.hg38.bed.gz"
    if repeat_gz.exists():
        extract_line_l1(repeat_gz, line_l1, args.force)
    else:
        print(f"WARNING: missing {repeat_gz}")

    snp_win = SUPPORT / f"common_snp_windows_pm{window_bp}bp.hg38.bed.gz"
    if not args.skip_snp_windows:
        # Prefer already-filtered meqtl pvar (smaller, MAF-filtered)
        cand = [
            PROJECT_ROOT / "meqtl-validation/01_cpg_meqtl_mapping/caudate/_m/genotypes/meqtl_AA.pvar",
            PROJECT_ROOT / paths["genotype"]["AA"]["pvar"],
        ]
        pvar = next((p for p in cand if Path(p).exists()), None)
        if pvar is None:
            print("WARNING: no pvar found for SNP proximity windows")
        else:
            build_snp_proximity_windows(Path(pvar), snp_win, maf_min, window_bp, args.force)

    assets = []
    for path, url, notes in [
        (blacklist, "Boyle-Lab Blacklist hg38 v2", "exclude_blacklist"),
        (segdup_txt, "UCSC genomicSuperDups.txt.gz", "source_table"),
        (segdup_bed, "converted from genomicSuperDups", "exclude_segmental_duplications"),
        (umap_bw, "Hoffman k24 Umap MultiTrack", "high_mappability_scoring"),
        (ct_bed, str(CT_SNP_DIR), "exclude_common_ct_snp_at_cpg"),
        (line_l1, "subset of repeat-masker-hg38.gz", "LINE_L1_enrichment"),
        (snp_win, "from meqtl_AA.pvar or TOPMed AA pvar", f"exclude_snp_proximal_cpgs_{window_bp}bp"),
        (repeat_gz, "existing project asset", "repeat_overlap"),
    ]:
        if not path.exists() or path.stat().st_size == 0:
            assets.append({
                "filename": path.name,
                "path": str(path),
                "url_or_source": url,
                "bytes": 0,
                "sha256": "",
                "download_utc": now,
                "status": "missing",
                "notes": notes,
            })
            continue
        assets.append({
            "filename": path.name,
            "path": str(path),
            "url_or_source": url,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "download_utc": now,
            "status": "ready",
            "notes": notes,
        })

    write_tsv(SUPPORT / "annotation_asset_manifest.tsv", assets)
    # also mirror under Phase 6 module
    phase6 = PROJECT_ROOT / "meqtl-validation/07_repeat_mappability_sensitivity/_m"
    phase6.mkdir(parents=True, exist_ok=True)
    write_tsv(phase6 / "annotation_asset_manifest.tsv", assets)
    print(f"Wrote manifest with {len(assets)} assets")


if __name__ == "__main__":
    main()
