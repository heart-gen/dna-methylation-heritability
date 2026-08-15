#!/usr/bin/env python3
"""
Convert staff-repo CpG phen matrices to TensorQTL phenotype BED for VMR CpGs.

Input phen format (Alexis repo):
  header: FID IID <cpg_pos> <cpg_pos> ...
  rows: BrNum  arrayID  beta_value ...

Output TensorQTL BED (bgzip + tabix):
  #chr start end phenotype_id sample1 sample2 ...
  Positions are 0-based half-open for a single CpG (pos-1, pos).
  Phenotype IDs: {chrom}_{pos} (1-based genomic position of the C).

Only CpGs overlapping the region VMR BED are retained. Samples are restricted
to the primary inclusion list from 00_preflight.py.
"""

from __future__ import annotations

import argparse
import gzip, sys
import subprocess
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import (  # noqa: E402
    cpg_matrix_relpath,
    load_paths,
    load_yaml,
    read_tsv,
    write_tsv,
)

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True, choices=["caudate", "dlpfc", "hippocampus"])
    p.add_argument("--chrom", required=True, help="Chromosome number without chr, e.g. 22")
    p.add_argument("--population", default="AA")
    p.add_argument(
        "--min-variance",
        type=float,
        default=0.0,
        help="Drop CpGs with sample variance <= this value among included samples",
    )
    return p.parse_args()


def load_vmrs(bed: Path, chrom: str, pred_path: Path | None = None) -> list[tuple[int, int, str]]:
    chrom_norm = str(chrom).replace("chr", "")
    pred_intervals: list[tuple[int, int, str]] = []
    if pred_path is not None and pred_path.exists():
        import csv
        with pred_path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                c = str(row.get("chrom", "")).replace("chr", "")
                if c != chrom_norm:
                    continue
                pred_intervals.append((int(row["start"]), int(row["end"]), str(row.get("task_id", ""))))
        pred_intervals.sort()

    def best_task_id(start: int, end: int) -> str | None:
        hits = []
        for ps, pe, tid in pred_intervals:
            if pe <= start:
                continue
            if ps >= end:
                break
            if start < pe and end > ps:
                overlap = min(end, pe) - max(start, ps)
                hits.append((overlap, tid))
        if not hits:
            return None
        hits.sort(reverse=True)
        return hits[0][1] or None

    intervals = []
    with bed.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[0].replace("chr", "") != chrom_norm:
                continue
            start, end = int(parts[1]), int(parts[2])
            tid = best_task_id(start, end)
            if tid:
                vid = tid
            elif len(parts) > 3 and parts[3]:
                vid = parts[3]
            else:
                vid = f"{chrom_norm}:{start}-{end}"
            intervals.append((start, end, vid))
    intervals.sort()
    return intervals


def overlaps_any(pos0: int, intervals: list[tuple[int, int, str]]) -> str | None:
    """pos0 is 0-based CpG start; intervals are half-open BED."""
    # linear scan is fine per chrom after filtering; use two-pointer for speed
    for start, end, vid in intervals:
        if pos0 >= end:
            continue
        if pos0 + 2 <= start:  # CpG spans 2 bp
            return None if start > pos0 + 2 else None
        if start <= pos0 < end or start < pos0 + 2 <= end or (pos0 <= start and pos0 + 2 >= end):
            return vid
        if pos0 + 2 < start:
            break
    # fallback simple
    for start, end, vid in intervals:
        if pos0 < end and pos0 + 2 > start:
            return vid
    return None


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def load_exclusion_intervals(path: Path, chrom: str) -> list[tuple[int, int]]:
    """Load one chromosome of a BED exclusion asset."""
    if not path.exists():
        raise SystemExit(f"Required exclusion asset is missing: {path}")
    wanted = f"chr{str(chrom).replace('chr', '')}"
    rows = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] == wanted:
                rows.append((int(parts[1]), int(parts[2])))
    rows.sort()
    return rows


def overlaps_interval(start: int, end: int, intervals: list[tuple[int, int]]) -> bool:
    for left, right in intervals:
        if left >= end:
            return False
        if right > start:
            return True
    return False


def main() -> None:
    args = parse_args()
    paths = load_paths()
    meqtl = load_yaml("meqtl_parameters.yml")
    project = Path(paths["project_root"])
    staff = Path(paths["cpg_methylation_root"])
    region = args.region
    chrom = str(args.chrom).replace("chr", "")

    preflight_dir = (
        project / "meqtl-validation" / "01_cpg_meqtl_mapping" / region / "_m" / "preflight"
    )
    if args.population != "AA":
        preflight_dir = preflight_dir / args.population
    preflight = preflight_dir / "sample_inclusion_primary.tsv"
    if not preflight.exists():
        raise SystemExit(f"Missing {preflight}; run 00_preflight.py first")
    keep = [r["brnum"] for r in read_tsv(preflight)]
    keep_set = set(keep)

    phen_path = staff / cpg_matrix_relpath(paths, region, chrom, args.population)
    if not phen_path.exists():
        raise SystemExit(f"Missing CpG matrix: {phen_path}")

    if args.population == "AA":
        pred_key = "local_predictability_summary_template"
    else:
        pred_key = "local_predictability_ea_template"
    vmr_key = "vmr_bed_aa_template" if args.population == "AA" else "vmr_bed_all_individuals_template"
    vmr_bed = project / paths.get(vmr_key, paths["vmr_bed_template"]).format(region=region)
    pred_path = project / paths[pred_key].format(region=region)
    intervals = load_vmrs(vmr_bed, chrom, pred_path=pred_path)
    if not intervals:
        print(f"No VMRs on chr{chrom}; writing empty outputs")
    
    prepared_root = (
        project / "meqtl-validation" / "01_cpg_meqtl_mapping" / region / "_m" / "prepared"
    )
    outdir = prepared_root if args.population == "AA" else prepared_root / args.population
    outdir.mkdir(parents=True, exist_ok=True)

    # Read phen: header positions, rows samples
    with phen_path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        site_cols = header[2:]  # genomic 1-based positions as strings
        sample_rows = []
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            fid = parts[0]
            if fid not in keep_set:
                continue
            sample_rows.append((fid, parts[2:]))

    # Order samples as in keep list for reproducibility
    by_fid = {fid: vals for fid, vals in sample_rows}
    ordered = [br for br in keep if br in by_fid]
    if not ordered:
        raise SystemExit("No overlapping samples between inclusion list and CpG matrix")

    blacklist = []
    if meqtl["cpg_qc"].get("exclude_blacklist", False):
        blacklist = load_exclusion_intervals(
            project / paths["support_files"]["blacklist_hg38"], chrom
        )
    ct_snps = []
    if meqtl["cpg_qc"].get("exclude_common_ct_snp_at_cpg", False):
        ct_snps = load_exclusion_intervals(
            project / paths["support_files"]["ct_snps_at_cpg"], chrom
        )

    # Select CpG columns overlapping VMRs after prespecified exclusions.
    selected = []  # (col_idx, pos1, vmr_id)
    n_in_vmr_before_exclusion = 0
    n_excluded_blacklist = 0
    n_excluded_ct_snp = 0
    for idx, pos_str in enumerate(site_cols):
        try:
            pos1 = int(float(pos_str))
        except ValueError:
            continue
        pos0 = pos1 - 1
        vid = None
        for start, end, v in intervals:
            if pos0 < end and pos0 + 2 > start:
                vid = v
                break
        if vid is None:
            continue
        n_in_vmr_before_exclusion += 1
        if blacklist and overlaps_interval(pos0, pos0 + 2, blacklist):
            n_excluded_blacklist += 1
            continue
        if ct_snps and overlaps_interval(pos0, pos0 + 2, ct_snps):
            n_excluded_ct_snp += 1
            continue
        selected.append((idx, pos1, vid))

    # Build matrix and filter variance
    mat = np.array(
        [[float(by_fid[br][idx]) if by_fid[br][idx] not in ("", "NA", "NaN") else np.nan
          for idx, _, _ in selected]
         for br in ordered],
        dtype=float,
    )
    keep_cpg_mask = []
    cpg_map_rows = []
    chrom_label = f"chr{chrom}"
    for j, (idx, pos1, vid) in enumerate(selected):
        col = mat[:, j]
        n_nonmiss = np.isfinite(col).sum()
        frac = n_nonmiss / len(ordered) if ordered else 0
        var = float(np.nanvar(col)) if n_nonmiss > 1 else 0.0
        mean_methylation = float(np.nanmean(col)) if n_nonmiss else np.nan
        min_frac = meqtl["cpg_qc"]["min_fraction_samples_passing_coverage"]
        # Phen matrices are prefiltered on depth in Alexis stats.rda pipeline
        # (cov >= min_coverage in >= min_frac samples). Non-missing fraction is
        # retained as a secondary completeness check; real depth is attached
        # later by 06_extract_cpg_coverage.R (mean_coverage, frac cov>=min).
        ok = frac >= min_frac and var > args.min_variance
        if meqtl["cpg_qc"]["require_nonzero_variance"] and var <= 0:
            ok = False
        keep_cpg_mask.append(ok)
        if ok:
            pid = f"{chrom_label}_{pos1}"
            cpg_map_rows.append({
                "phenotype_id": pid,
                "chrom": chrom_label,
                "pos_1based": pos1,
                "start_0based": pos1 - 1,
                "end_0based": pos1,
                "vmr_id": vid,
                "n_nonmissing": int(n_nonmiss),
                "fraction_nonmissing": round(frac, 4),
                "variance": var,
                "mean_methylation": mean_methylation,
                "coverage_note": "run_06_extract_cpg_coverage_for_depth",
            })

    keep_idx = [i for i, ok in enumerate(keep_cpg_mask) if ok]
    bed_path = outdir / f"cpg_phenotypes.{chrom_label}.bed"
    with bed_path.open("w") as handle:
        handle.write("\t".join(["#chr", "start", "end", "phenotype_id"] + ordered) + "\n")
        for j in keep_idx:
            idx, pos1, vid = selected[j]
            pid = f"{chrom_label}_{pos1}"
            start0 = pos1 - 1
            vals = []
            for br in ordered:
                raw = by_fid[br][idx]
                vals.append("NA" if raw in ("", "NA", "NaN") else raw)
            handle.write("\t".join([chrom_label, str(start0), str(pos1), pid] + vals) + "\n")

    # bgzip + tabix if available
    gz_path = Path(str(bed_path) + ".gz")
    try:
        subprocess.check_call(["bgzip", "-f", str(bed_path)])
        subprocess.check_call(["tabix", "-p", "bed", str(gz_path)])
        bed_out = gz_path
    except (FileNotFoundError, subprocess.CalledProcessError):
        # leave uncompressed bed
        bed_out = bed_path
        print("WARNING: bgzip/tabix unavailable; left uncompressed BED")

    write_tsv(
        outdir / f"cpg_vmr_map.{chrom_label}.tsv",
        cpg_map_rows,
        [
            "phenotype_id", "chrom", "pos_1based", "start_0based",
            "end_0based", "vmr_id", "n_nonmissing", "fraction_nonmissing",
            "variance", "mean_methylation", "coverage_note",
        ],
    )
    write_tsv(
        outdir / f"prep_summary.{chrom_label}.tsv",
        [{
            "region": region,
            "chrom": chrom_label,
            "n_samples": len(ordered),
            "n_sites_in_phen": len(site_cols),
            "n_sites_in_vmrs": len(selected),
            "n_sites_in_vmrs_before_exclusion": n_in_vmr_before_exclusion,
            "n_excluded_blacklist": n_excluded_blacklist,
            "n_excluded_common_ct_snp": n_excluded_ct_snp,
            "n_sites_retained": len(keep_idx),
            "bed_path": str(bed_out),
        }],
    )
    write_tsv(outdir / f"sample_order.{chrom_label}.tsv", [{"brnum": b} for b in ordered], ["brnum"])
    print(f"Wrote {bed_out} with {len(keep_idx)} CpGs for {len(ordered)} samples")


if __name__ == "__main__":
    main()
