#!/usr/bin/env python3
"""Reaggregate locked M3a-tested CpGs into all-individual VMR coordinates.

This is a coordinate-definition sensitivity. It never treats an untested CpG as
nonsignificant: the denominator contains only CpGs returned by locked M3a.
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import ANALYSIS_SCHEMA_VERSION, canonical_vmr_id, load_paths, write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", required=True)
    parser.add_argument("--lead", required=True)
    parser.add_argument("--vmr-bed", required=True)
    parser.add_argument("--prepared-map-glob", required=True)
    parser.add_argument("--legacy", required=True)
    parser.add_argument("--technical", required=True)
    parser.add_argument("--annotation", default="")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args()


def _chrom(value: object) -> str:
    return str(value).replace("chr", "")


def _read_vmrs(path: Path) -> pd.DataFrame:
    vmrs = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    vmrs["chrom"] = vmrs["chrom"].map(_chrom)
    vmrs[["start", "end"]] = vmrs[["start", "end"]].astype(int)
    vmrs = vmrs.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    vmrs["vmr_coord_id"] = [
        canonical_vmr_id(chrom, start, end)
        for chrom, start, end in zip(vmrs["chrom"], vmrs["start"], vmrs["end"])
    ]
    if vmrs["vmr_coord_id"].duplicated().any():
        raise SystemExit("All-individual VMR BED contains duplicate coordinates")
    for chrom, group in vmrs.groupby("chrom", sort=False):
        ordered = group.sort_values("start")
        if (ordered["start"].to_numpy()[1:] <= ordered["end"].to_numpy()[:-1]).any():
            raise SystemExit(f"Overlapping all-individual VMRs on chromosome {chrom}")
    return vmrs


def assign_coordinates(points: pd.DataFrame, vmrs: pd.DataFrame) -> pd.Series:
    """Assign 1-based CpG positions to inclusive project VMR coordinates."""
    assigned = pd.Series(pd.NA, index=points.index, dtype="object")
    for chrom, point_group in points.groupby("chrom", sort=False):
        intervals = vmrs[vmrs["chrom"].eq(chrom)].sort_values("start")
        if intervals.empty:
            continue
        starts = intervals["start"].to_numpy(dtype=int)
        ends = intervals["end"].to_numpy(dtype=int)
        ids = intervals["vmr_coord_id"].to_numpy(dtype=object)
        positions = point_group["pos_1based"].to_numpy(dtype=int)
        loc = np.searchsorted(starts, positions, side="right") - 1
        valid = loc >= 0
        clipped = np.clip(loc, 0, len(starts) - 1)
        valid &= positions <= ends[clipped]
        assigned.loc[point_group.index[valid]] = ids[clipped[valid]]
    return assigned


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _aggregate(tested: pd.DataFrame, region: str, fdr: float) -> pd.DataFrame:
    tested = tested.copy()
    tested["sig_meqtl"] = _numeric(tested, "qval").le(fdr)
    tested["_p"] = _numeric(tested, "pval_beta")
    if tested["_p"].isna().all():
        tested["_p"] = _numeric(tested, "pval_nominal")
    tested["_slope"] = _numeric(tested, "slope").abs()
    tested["_distance"] = _numeric(tested, "start_distance").abs()
    tested["_cov"] = _numeric(tested, "mean_coverage")
    tested["_var"] = _numeric(tested, "variance")
    tested["_mean"] = _numeric(tested, "mean_methylation")
    tested["_num_var"] = _numeric(tested, "num_var")
    tested["_sig_variant"] = tested.get("variant_id", pd.Series(pd.NA, index=tested.index)).where(
        tested["sig_meqtl"]
    )
    result = tested.groupby("vmr_coord_id", as_index=False).agg(
        n_tested_cpgs=("phenotype_id", "size"),
        n_cpgs_with_sig_meqtl=("sig_meqtl", "sum"),
        min_cpg_pvalue=("_p", "min"),
        strongest_abs_beta=("_slope", "max"),
        median_cpg_lead_snp_distance=("_distance", "median"),
        average_cpg_coverage=("_cov", "mean"),
        mean_cpg_variance=("_var", "mean"),
        vmr_mean_methylation=("_mean", "mean"),
        mean_num_tested_snps_per_cpg=("_num_var", "mean"),
        n_distinct_significant_lead_snps=("_sig_variant", "nunique"),
    )
    result["n_cpgs_with_sig_meqtl"] = result["n_cpgs_with_sig_meqtl"].astype(int)
    result["proportion_cpgs_with_sig_meqtl"] = (
        result["n_cpgs_with_sig_meqtl"] / result["n_tested_cpgs"]
    )
    result["region"] = region
    return result


def _annotation(path: Path) -> pd.DataFrame:
    if not path.is_file():
        return pd.DataFrame(columns=["vmr_coord_id", "genomic_annotation"])
    frame = pd.read_csv(path, sep="\t")
    required = {"seqnames", "start", "end"}
    if not required.issubset(frame.columns):
        return pd.DataFrame(columns=["vmr_coord_id", "genomic_annotation"])
    frame["vmr_coord_id"] = [
        canonical_vmr_id(chrom, start, end)
        for chrom, start, end in zip(frame["seqnames"], frame["start"], frame["end"])
    ]
    label = "annot.type" if "annot.type" in frame.columns else None
    if label is None:
        return pd.DataFrame(columns=["vmr_coord_id", "genomic_annotation"])
    return frame.groupby("vmr_coord_id", as_index=False).agg(
        genomic_annotation=(label, lambda values: ";".join(sorted(set(values.dropna().astype(str)))))
    )


def main() -> None:
    args = parse_args()
    region = args.region.lower()
    lead_path = Path(args.lead)
    vmr_path = Path(args.vmr_bed)
    map_paths = [Path(path) for path in sorted(glob.glob(args.prepared_map_glob))]
    for path in [lead_path, vmr_path, Path(args.legacy), Path(args.technical)]:
        if not path.is_file():
            raise SystemExit(f"Required input is missing: {path}")
    if not map_paths:
        raise SystemExit(f"No prepared maps match {args.prepared_map_glob}")

    vmrs = _read_vmrs(vmr_path)
    metadata = pd.concat([pd.read_csv(path, sep="\t") for path in map_paths], ignore_index=True)
    metadata = metadata.drop_duplicates("phenotype_id", keep="first")
    metadata["chrom"] = metadata["chrom"].map(_chrom)
    metadata["pos_1based"] = pd.to_numeric(metadata["pos_1based"], errors="raise").astype(int)
    metadata["aligned_vmr_coord_id"] = assign_coordinates(metadata, vmrs)

    lead = pd.read_csv(lead_path, sep="\t")
    if lead["phenotype_id"].duplicated().any():
        raise SystemExit("Locked lead table is not one row per tested CpG")
    tested = lead.merge(
        metadata.drop(columns=["vmr_id"], errors="ignore"),
        on="phenotype_id", how="left", validate="one_to_one", suffixes=("", "_prepared"),
    )
    tested["vmr_coord_id"] = tested["aligned_vmr_coord_id"]
    n_tested_total = len(tested)
    tested = tested.dropna(subset=["vmr_coord_id"]).copy()
    aggregated = _aggregate(tested, region, float(args.fdr))

    prepared_counts = (
        metadata.dropna(subset=["aligned_vmr_coord_id"])
        .groupby("aligned_vmr_coord_id").size().rename("n_prepared_m3a_cpgs").reset_index()
        .rename(columns={"aligned_vmr_coord_id": "vmr_coord_id"})
    )
    burden = vmrs.merge(prepared_counts, on="vmr_coord_id", how="left")
    burden = burden.merge(aggregated, on="vmr_coord_id", how="left")
    for column in ["n_prepared_m3a_cpgs", "n_tested_cpgs", "n_cpgs_with_sig_meqtl", "n_distinct_significant_lead_snps"]:
        burden[column] = burden[column].fillna(0).astype(int)
    burden["region"] = burden["region"].fillna(region)
    burden["n_prepared_but_untested_cpgs"] = burden["n_prepared_m3a_cpgs"] - burden["n_tested_cpgs"]
    burden["vmr_id"] = burden["vmr_coord_id"]
    burden["length"] = burden["end"] - burden["start"]
    burden["cpg_density"] = burden["n_tested_cpgs"] / burden["length"].replace(0, np.nan)

    legacy = pd.read_csv(args.legacy, sep="\t")
    legacy["vmr_coord_id"] = [
        canonical_vmr_id(chrom, start, end)
        for chrom, start, end in zip(legacy["chrom"], legacy["start"], legacy["end"])
    ]
    legacy = legacy.drop_duplicates("vmr_coord_id", keep="first")
    burden = burden.merge(
        legacy[["vmr_coord_id", "h2_unscaled", "num_snps"]].rename(
            columns={"h2_unscaled": "local_predictability"}
        ), on="vmr_coord_id", how="left", validate="one_to_one",
    )
    burden["local_common_snp_density"] = burden["num_snps"] / 1_000_001.0

    technical = pd.read_csv(args.technical, sep="\t")
    technical["vmr_coord_id"] = [
        canonical_vmr_id(chrom, start, end)
        for chrom, start, end in zip(technical["chrom"], technical["start"], technical["end"])
    ]
    tech_columns = [column for column in [
        "umap_k24_mean", "high_mappability", "line_l1_frac", "segdup_frac",
        "overlaps_segdup", "blacklist_frac", "overlaps_blacklist",
        "snp_prox_frac", "overlaps_snp_prox",
    ] if column in technical.columns]
    technical = technical.drop_duplicates("vmr_coord_id", keep="first")
    burden = burden.merge(
        technical[["vmr_coord_id"] + tech_columns], on="vmr_coord_id", how="left", validate="one_to_one"
    )
    burden["tech_join_source"] = np.where(
        burden[tech_columns[0]].notna() if tech_columns else False, "exact_coord", pd.NA
    )
    burden = burden.merge(_annotation(Path(args.annotation)), on="vmr_coord_id", how="left", validate="one_to_one")
    burden["analysis_schema_version"] = ANALYSIS_SCHEMA_VERSION
    burden["coordinate_sensitivity"] = "all_individual_vmr_locked_m3a_tested_cpgs_only"

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    burden_path = output / "vmr_meqtl_burden.tsv.gz"
    burden.to_csv(burden_path, sep="\t", index=False, compression="gzip")
    write_tsv(output / "coordinate_alignment_qc.tsv", [{
        "region": region,
        "analysis_schema_version": ANALYSIS_SCHEMA_VERSION,
        "n_all_individual_vmrs": int(len(vmrs)),
        "n_locked_m3a_tested_cpgs": int(n_tested_total),
        "n_tested_cpgs_assigned": int(len(tested)),
        "n_tested_cpgs_unassigned": int(n_tested_total - len(tested)),
        "fraction_tested_cpgs_assigned": float(len(tested) / n_tested_total),
        "n_vmrs_with_tested_cpgs": int(burden["n_tested_cpgs"].gt(0).sum()),
        "n_vmrs_with_legacy": int(burden["local_predictability"].notna().sum()),
        "n_vmrs_with_technical_annotations": int(burden[tech_columns[0]].notna().sum()) if tech_columns else 0,
        "coordinate_convention": "project VMR start/end; CpG pos_1based assigned inclusively",
        "denominator_note": "only CpGs returned by locked M3a are counted as tested",
        "output": str(burden_path),
    }])
    print(f"Wrote coordinate-aligned burden to {burden_path}")


if __name__ == "__main__":
    main()
