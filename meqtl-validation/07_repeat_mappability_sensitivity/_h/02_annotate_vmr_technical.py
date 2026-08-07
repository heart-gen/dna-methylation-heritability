#!/usr/bin/env python3
"""
Annotate VMRs with technical assets for Phase 6 sensitivity analyses.

Outputs per region:
  vmr_technical_annotations.tsv
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import PROJECT_ROOT, load_paths, load_yaml, write_tsv  # noqa: E402

SUPPORT = PROJECT_ROOT / "inputs" / "supportfiles" / "_m"
OUTDIR = PROJECT_ROOT / "meqtl-validation/07_repeat_mappability_sensitivity/_m"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions", nargs="+", default=["caudate", "dlpfc", "hippocampus"])
    return p.parse_args()


def load_vmr_bed(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0] if parts[0].startswith("chr") else f"chr{parts[0]}"
            start, end = int(parts[1]), int(parts[2])
            rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "length": end - start,
                "interval_id": f"{chrom.replace('chr', '')}:{start}-{end}",
            })
    return pd.DataFrame(rows)


def attach_task_ids(vmrs: pd.DataFrame, pred_path: Path) -> pd.DataFrame:
    """Mirror Phase 1 VMR id choice: use elastic-net task_id when interval matches."""
    out = vmrs.copy()
    if not pred_path.exists():
        out["task_id"] = pd.NA
        out["vmr_id"] = out["interval_id"]
        return out
    p = pd.read_csv(pred_path, sep="\t")
    if not {"chrom", "start", "end", "task_id"}.issubset(p.columns):
        out["task_id"] = pd.NA
        out["vmr_id"] = out["interval_id"]
        return out
    p = p.copy()
    p["chrom"] = p["chrom"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")
    p["start"] = p["start"].astype(int)
    p["end"] = p["end"].astype(int)
    p["task_id"] = p["task_id"].map(lambda x: str(int(x)) if pd.notna(x) else pd.NA)
    out = out.merge(p[["chrom", "start", "end", "task_id"]], on=["chrom", "start", "end"], how="left")
    out["vmr_id"] = out["task_id"].where(out["task_id"].notna(), out["interval_id"])
    return out


def bedtools_coverage_fraction(vmrs: pd.DataFrame, feature_bed: Path) -> pd.Series:
    """Return overlap fraction aligned to vmrs.index using unique bed names."""
    n = len(vmrs)
    if not feature_bed.exists() or feature_bed.stat().st_size == 0:
        return pd.Series([pd.NA] * n, index=vmrs.index, dtype="float")
    with tempfile.TemporaryDirectory(prefix="vmr_ann_") as tmp:
        tmp_path = Path(tmp)
        work = vmrs.copy()
        work["name"] = [f"v{i}" for i in range(n)]
        vmr_bed = tmp_path / "vmr.bed"
        work[["chrom", "start", "end", "name"]].to_csv(vmr_bed, sep="\t", header=False, index=False)
        vmr_sorted = tmp_path / "vmr.sorted.bed"
        feat_sorted = tmp_path / "feat.sorted.bed"
        with vmr_sorted.open("w") as hout:
            subprocess.check_call(["bedtools", "sort", "-i", str(vmr_bed)], stdout=hout)
        if str(feature_bed).endswith(".gz"):
            feat_raw = tmp_path / "feat.bed"
            with gzip.open(feature_bed, "rt") as hin, feat_raw.open("w") as hout:
                for line in hin:
                    if line.strip() and not line.startswith("#"):
                        hout.write(line)
            feat_in = feat_raw
        else:
            feat_in = feature_bed
        with feat_sorted.open("w") as hout:
            subprocess.check_call(["bedtools", "sort", "-i", str(feat_in)], stdout=hout)
        out = tmp_path / "cov.tsv"
        with out.open("w") as hout:
            subprocess.check_call(
                ["bedtools", "coverage", "-a", str(vmr_sorted), "-b", str(feat_sorted)],
                stdout=hout,
            )
        cov = pd.read_csv(
            out, sep="\t", header=None,
            names=["chrom", "start", "end", "name", "n_features", "nbases", "length", "frac"],
            dtype={"name": str},
        )
        frac_map = cov.set_index("name")["frac"]
        return work["name"].map(frac_map)


def mean_umap(vmrs: pd.DataFrame, bw: Path) -> pd.Series:
    tool = shutil.which("bigWigAverageOverBed")
    if tool is None or not bw.exists():
        print("WARNING: bigWigAverageOverBed or Umap bigWig missing; mappability left NA")
        return pd.Series([pd.NA] * len(vmrs), index=vmrs.index, dtype="float")
    with tempfile.TemporaryDirectory(prefix="umap_") as tmp:
        tmp_path = Path(tmp)
        bed = tmp_path / "vmr.bed"
        out = tmp_path / "avg.tab"
        work = vmrs.copy()
        work["name"] = [f"v{i}" for i in range(len(work))]
        work[["chrom", "start", "end", "name"]].to_csv(bed, sep="\t", header=False, index=False)
        subprocess.check_call([tool, str(bw), str(bed), str(out)])
        avg = pd.read_csv(
            out, sep="\t", header=None,
            names=["name", "size", "covered", "sum", "mean0", "mean"],
            dtype={"name": str},
        )
        return work["name"].map(avg.set_index("name")["mean"])


def annotate_region(region: str, paths: dict, thresholds: dict) -> Path:
    vmr_bed = PROJECT_ROOT / paths["vmr_bed_template"].format(region=region)
    pred = PROJECT_ROOT / paths["local_predictability_summary_template"].format(region=region)
    vmrs = attach_task_ids(load_vmr_bed(vmr_bed), pred)

    blacklist = SUPPORT / "hg38-blacklist.v2.bed.gz"
    segdup = SUPPORT / "genomicSuperDups.hg38.bed.gz"
    line_l1 = SUPPORT / "repeat-masker.LINE_L1.hg38.bed.gz"
    umap_bw = SUPPORT / "k24.Umap.MultiTrackMappability.bw"

    out = vmrs.copy()
    out["blacklist_frac"] = bedtools_coverage_fraction(vmrs, blacklist).astype(float)
    out["segdup_frac"] = bedtools_coverage_fraction(vmrs, segdup).astype(float)
    out["line_l1_frac"] = bedtools_coverage_fraction(vmrs, line_l1).astype(float)
    out["umap_k24_mean"] = mean_umap(vmrs, umap_bw).astype(float)
    out["high_mappability"] = (
        out["umap_k24_mean"] >= float(thresholds["repeat_sensitivity"]["high_mappability_min"])
    ).astype(int)
    out["overlaps_blacklist"] = (out["blacklist_frac"].fillna(0) > 0).astype(int)
    out["overlaps_segdup"] = (out["segdup_frac"].fillna(0) > 0).astype(int)
    snp_win = SUPPORT / "common_snp_windows_pm150bp.hg38.bed.gz"
    if snp_win.exists():
        out["snp_prox_frac"] = bedtools_coverage_fraction(vmrs, snp_win).astype(float)
        out["overlaps_snp_prox"] = (out["snp_prox_frac"].fillna(0) > 0).astype(int)
    out["region"] = region

    dest = OUTDIR / region
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / "vmr_technical_annotations.tsv"
    cols = [
        "region", "chrom", "start", "end", "length", "vmr_id", "interval_id", "task_id",
        "blacklist_frac", "segdup_frac", "line_l1_frac", "umap_k24_mean",
        "high_mappability", "overlaps_blacklist", "overlaps_segdup",
    ]
    if "snp_prox_frac" in out.columns:
        cols += ["snp_prox_frac", "overlaps_snp_prox"]
    out[cols].to_csv(path, sep="\t", index=False)
    print(
        f"{region}: wrote {path} ({len(out):,} VMRs; "
        f"high_mappability={int(out['high_mappability'].sum())}; "
        f"blacklist_NA={int(out['blacklist_frac'].isna().sum())})"
    )
    return path


def main() -> None:
    args = parse_args()
    paths = load_paths()
    thresholds = load_yaml("analysis_thresholds.yml")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for region in args.regions:
        p = annotate_region(region, paths, thresholds)
        rows.append({"region": region, "path": str(p), "status": "ready"})
    write_tsv(OUTDIR / "vmr_technical_annotation_index.tsv", rows)


if __name__ == "__main__":
    main()
