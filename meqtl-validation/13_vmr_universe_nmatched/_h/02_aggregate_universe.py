#!/usr/bin/env python3
"""
Aggregate the per-chromosome universe decomposition into a verdict.

Reads the 66 outputs of 01_universe_decomposition.R and answers: is caudate's
larger testable VMR/CpG universe a sample-size artifact, or a property of the
data?

The comparison ladder, in order of what each rung controls for:

  1. WITHIN-REGION, N ONLY.  caudate at N=111 (30 module-10 replicate donor
     lists) vs caudate at full N. Same tissue, same libraries, same CpGs --
     the only thing that changes is how many donors enter the 80%-of-donors
     coverage filter. This isolates sample size exactly, with nothing else to
     adjust for.

  2. BETWEEN-REGION, N MATCHED.  caudate at N=111 vs DLPFC at N=111. Any gap
     surviving here is by construction not sample size.

  3. WHY, IF NOT N.  Per-donor coverage by region, restricted to the 92 donors
     assayed in all three regions. Pairing on donor removes donor-level
     variation entirely, so a residual difference is a library/sequencing
     property of the region's data rather than anything about the cohort.

Rung 3 is the one that carries interpretive weight, and it separates two
technical explanations that are easy to conflate:

  DEPTH       caudate libraries sequenced deeper -> more CpGs clear cov>=5.
  UNIFORMITY  caudate libraries more consistent between donors -> more CpGs
              clear cov>=5 in >=80% of donors, even at lower mean depth.

The filter is a statement about the whole donor distribution at a CpG, so it is
broken by the low tail rather than by the mean. A shallower but tighter region
keeps more CpGs. Report both; do not read mean depth alone.

Usage:
    conda activate /projects/p32505/opt/envs/genomics
    python3 ../_h/02_aggregate_universe.py
"""

from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import ANALYSIS_SCHEMA_VERSION, write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
MODULE = PROJECT / "meqtl-validation/13_vmr_universe_nmatched"
INDIR = MODULE / "_m"
REGIONS = ("caudate", "dlpfc", "hippocampus")
AUTOSOMES = [str(c) for c in range(1, 23)]
SHARED3 = (
    PROJECT
    / "meqtl-validation/04_cross_region_sharing/_m/caudate_downsample/shared3_donors.txt"
)

# Full-N donor counts, from module 10's design summary.
DESIGN_N = {"caudate": 153, "dlpfc": 111, "hippocampus": 116}
# A between-region gap this small is not worth interpreting as a real
# difference in testable universe.
NEGLIGIBLE_RATIO = 0.02


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--indir", default=str(INDIR))
    p.add_argument("--outdir", default=str(INDIR))
    p.add_argument(
        "--require-complete",
        action="store_true",
        help="Fail if any region x autosome output is missing, rather than reporting partial results.",
    )
    return p.parse_args()


def load_kind(indir: Path, kind: str) -> pd.DataFrame:
    """Concatenate all per-region, per-chromosome files of one output kind."""
    frames = []
    missing = []
    for region in REGIONS:
        for chrom in AUTOSOMES:
            path = indir / f"{kind}.{region}.chr{chrom}.tsv"
            if not path.exists():
                missing.append(f"{region}:chr{chrom}")
                continue
            frames.append(pd.read_csv(path, sep="\t", dtype={"chr": str}))
    if not frames:
        raise SystemExit(f"No {kind} outputs found under {indir}")
    df = pd.concat(frames, ignore_index=True)
    df.attrs["missing"] = missing
    return df


def genome_totals(passes: pd.DataFrame) -> pd.DataFrame:
    """Sum CpG pass counts over autosomes, per region and per replicate.

    A replicate's genome-wide universe is the sum of its own per-chromosome
    counts, so replicates are summed within replicate id -- never averaged
    per chromosome and then summed, which would discard between-replicate
    variance and understate the spread.
    """
    full = (
        passes[passes["design"] == "full_n"]
        .groupby("region", as_index=False)
        .agg(n_sites_pass=("n_sites_pass", "sum"), n_sites_pre=("n_sites_pre", "sum"),
             n_chr=("chr", "nunique"), n_donors=("n_donors", "max"))
    )
    full["design"] = "full_n"
    full["replicate"] = pd.NA

    matched = passes[passes["design"] == "n_matched"]
    if matched.empty:
        return full

    per_rep = (
        matched.groupby(["region", "replicate"], as_index=False)
        .agg(n_sites_pass=("n_sites_pass", "sum"), n_sites_pre=("n_sites_pre", "sum"),
             n_chr=("chr", "nunique"), n_donors=("n_donors", "max"))
    )
    per_rep["design"] = "n_matched"
    return pd.concat([full, per_rep], ignore_index=True)


def summarize_rungs(totals: pd.DataFrame, coverage: pd.DataFrame) -> tuple[list[dict], dict]:
    rows: list[dict] = []
    verdict: dict = {}

    full = totals[totals["design"] == "full_n"].set_index("region")["n_sites_pass"].to_dict()
    matched = totals[(totals["design"] == "n_matched") & (totals["region"] == "caudate")]
    reps = sorted(matched["n_sites_pass"].tolist())

    # ---- rung 1: within caudate, sample size only ---------------------------
    if reps:
        mean_matched = statistics.fmean(reps)
        ratio_n = mean_matched / full["caudate"]
        rows.append(
            {
                "comparison": "caudate_n_matched_vs_caudate_full_n",
                "controls_for": "everything except donor count",
                "isolates": "sample size",
                "value_a": round(mean_matched, 1),
                "value_b": full["caudate"],
                "ratio": round(ratio_n, 4),
                "n_replicates": len(reps),
                "replicate_min": reps[0],
                "replicate_max": reps[-1],
                "full_n_within_replicate_range": bool(reps[0] <= full["caudate"] <= reps[-1]),
            }
        )
        verdict["sample_size_effect_ratio"] = round(ratio_n, 4)
        verdict["sample_size_explains_universe"] = bool(abs(ratio_n - 1.0) > NEGLIGIBLE_RATIO)
    else:
        verdict["sample_size_effect_ratio"] = None
        verdict["sample_size_explains_universe"] = None
        mean_matched = None

    # ---- rung 2: between region, N matched ----------------------------------
    for comparator in ("dlpfc", "hippocampus"):
        if comparator not in full:
            continue
        # Compare caudate's N-matched universe against the comparator's own
        # full-N universe: the comparator is already at (or near) N=111, so
        # this is the N-matched contrast.
        a = mean_matched if mean_matched is not None else full["caudate"]
        ratio = a / full[comparator]
        rows.append(
            {
                "comparison": f"caudate_n_matched_vs_{comparator}",
                "controls_for": "donor count",
                "isolates": "region / library",
                "value_a": round(a, 1),
                "value_b": full[comparator],
                "ratio": round(ratio, 4),
                "n_replicates": len(reps),
                "replicate_min": reps[0] if reps else "",
                "replicate_max": reps[-1] if reps else "",
                "full_n_within_replicate_range": "",
            }
        )
        verdict[f"n_matched_ratio_vs_{comparator}"] = round(ratio, 4)

    # ---- rung 3: per-donor coverage, paired on shared donors ----------------
    cov_rows, cov_verdict = coverage_contrast(coverage)
    rows.extend(cov_rows)
    verdict.update(cov_verdict)

    return rows, verdict


def coverage_contrast(coverage: pd.DataFrame) -> tuple[list[dict], dict]:
    """Per-donor mean coverage by region, unpaired and paired on shared donors."""
    rows: list[dict] = []
    verdict: dict = {}

    # Average each donor's coverage over autosomes, weighting chromosomes
    # equally -- we want a per-donor library property, not a per-base one.
    per_donor = (
        coverage.groupby(["region", "brnum"], as_index=False)
        .agg(mean_cov=("mean_cov", "mean"), frac_ge_min=("frac_ge_min", "mean"))
    )

    # Depth and uniformity are recorded separately and on purpose. The filter is
    # "at least 80% of donors have cov >= 5 here", which is a statement about
    # the whole donor distribution at a CpG, not about mean depth. A region can
    # be shallower on average and still keep more CpGs if its donors are
    # consistent, because the rule is broken by the low tail, not by the mean.
    for region in REGIONS:
        sub = per_donor[per_donor["region"] == region]
        if sub.empty:
            continue
        rows.append(
            {
                "comparison": f"per_donor_coverage_{region}",
                "controls_for": "-",
                "isolates": "library depth and uniformity",
                "value_a": round(sub["mean_cov"].mean(), 3),
                "value_b": round(sub["frac_ge_min"].mean(), 4),
                "ratio": "",
                "n_replicates": len(sub),
                "replicate_min": round(sub["mean_cov"].min(), 3),
                "replicate_max": round(sub["mean_cov"].max(), 3),
                "sd_mean_cov": round(sub["mean_cov"].std(), 3),
                "sd_frac_ge_min": round(sub["frac_ge_min"].std(), 4),
                "min_frac_ge_min": round(sub["frac_ge_min"].min(), 4),
                "full_n_within_replicate_range": "",
            }
        )
        verdict[f"mean_coverage_{region}"] = round(sub["mean_cov"].mean(), 3)
        verdict[f"sd_coverage_{region}"] = round(sub["mean_cov"].std(), 3)
        verdict[f"mean_frac_adequate_{region}"] = round(sub["frac_ge_min"].mean(), 4)
        verdict[f"sd_frac_adequate_{region}"] = round(sub["frac_ge_min"].std(), 4)

    # Paired on the 92 donors assayed in all three regions. Pairing removes
    # donor-level variation, so a residual gap is a property of the region's
    # libraries and not of which people happen to be in each cohort.
    if SHARED3.exists():
        shared = {line.strip() for line in SHARED3.read_text().splitlines() if line.strip()}
        wide = (
            per_donor[per_donor["brnum"].isin(shared)]
            .pivot(index="brnum", columns="region", values="mean_cov")
            .dropna()
        )
        verdict["n_shared_donors_paired"] = int(len(wide))
        for metric, label in (("mean_cov", "depth"), ("frac_ge_min", "adequacy")):
            wide_m = (
                per_donor[per_donor["brnum"].isin(shared)]
                .pivot(index="brnum", columns="region", values=metric)
                .dropna()
            )
            for comparator in ("dlpfc", "hippocampus"):
                if "caudate" not in wide_m.columns or comparator not in wide_m.columns:
                    continue
                delta = wide_m["caudate"] - wide_m[comparator]
                n_higher = int((delta > 0).sum())
                rows.append(
                    {
                        "comparison": f"paired_{label}_caudate_minus_{comparator}",
                        "controls_for": "donor identity",
                        "isolates": f"region library {label}",
                        "value_a": round(delta.mean(), 4),
                        "value_b": round(delta.median(), 4),
                        "ratio": round(wide_m["caudate"].mean() / wide_m[comparator].mean(), 4),
                        "n_replicates": len(wide_m),
                        "replicate_min": round(delta.min(), 4),
                        "replicate_max": round(delta.max(), 4),
                        "full_n_within_replicate_range": (
                            f"{n_higher}/{len(wide_m)} donors higher in caudate"
                        ),
                    }
                )
                verdict[f"paired_{label}_delta_vs_{comparator}"] = round(delta.mean(), 4)
                verdict[f"paired_{label}_frac_caudate_higher_vs_{comparator}"] = round(
                    n_higher / len(wide_m), 4
                )
    else:
        verdict["n_shared_donors_paired"] = 0

    return rows, verdict


def summarize_sd(sd: pd.DataFrame) -> list[dict]:
    """Mean residual-sd 99th percentile per region and design.

    02c.write_top_cpg.R thresholds at a within-region quantile, so this is the
    absolute variability bar each region's "top 1%" actually corresponds to.
    """
    out = (
        sd.groupby(["region", "design"], as_index=False)
        .agg(mean_sd_cutoff=("sd_cutoff", "mean"), n_chr=("chr", "nunique"),
             n_rows=("sd_cutoff", "size"))
        .sort_values(["region", "design"])
    )
    out["mean_sd_cutoff"] = out["mean_sd_cutoff"].round(6)
    return out.to_dict("records")


def build_verdict(verdict: dict, totals: pd.DataFrame) -> list[dict]:
    ss = verdict.get("sample_size_explains_universe")
    ratios = [
        verdict.get("n_matched_ratio_vs_dlpfc"),
        verdict.get("n_matched_ratio_vs_hippocampus"),
    ]
    ratios = [r for r in ratios if r is not None]
    surplus_survives = bool(ratios) and all(r > 1 + NEGLIGIBLE_RATIO for r in ratios)

    # Which technical property explains a surviving surplus: depth or uniformity?
    # These are separate, and the data can point at either. Read the paired
    # (donor-controlled) contrasts rather than the raw regional means.
    depth_deltas = [
        verdict.get("paired_depth_delta_vs_dlpfc"),
        verdict.get("paired_depth_delta_vs_hippocampus"),
    ]
    depth_deltas = [d for d in depth_deltas if d is not None]
    caudate_deeper = bool(depth_deltas) and all(d > 0 for d in depth_deltas)

    sd_c = verdict.get("sd_frac_adequate_caudate")
    sd_others = [
        verdict.get("sd_frac_adequate_dlpfc"),
        verdict.get("sd_frac_adequate_hippocampus"),
    ]
    sd_others = [s for s in sd_others if s is not None]
    caudate_more_uniform = (
        sd_c is not None and bool(sd_others) and all(sd_c < s for s in sd_others)
    )

    if ss is None:
        call = "incomplete"
        reading = "N-matched replicates absent; rerun step_1 for caudate."
    elif ss:
        call = "sample_size_artifact"
        reading = (
            "Caudate's universe shrinks materially when donors are N-matched, so the "
            "larger testable universe is a sample-size artifact and module 10's "
            "common-universe restriction is the correct comparison."
        )
    elif surplus_survives and caudate_more_uniform and not caudate_deeper:
        call = "coverage_uniformity_not_sample_size"
        reading = (
            "N-matching leaves caudate's universe essentially unchanged, so sample size "
            "explains none of it. Caudate is not sequenced deeper -- paired on shared "
            "donors it is shallower than both comparators -- but its per-donor coverage is "
            "markedly more uniform. The filter keeps a CpG when >=80% of donors clear "
            "cov>=5, which is broken by the low tail of the donor distribution rather than "
            "by mean depth, so the more consistent region retains more CpGs. The caudate "
            "universe surplus is a library-uniformity artifact: technical, not biological "
            "and not cohort size."
        )
    elif surplus_survives and caudate_deeper:
        call = "coverage_depth_not_sample_size"
        reading = (
            "N-matching leaves caudate's universe essentially unchanged, and caudate is "
            "sequenced deeper than both comparators on shared donors. The surplus is a "
            "sequencing-depth artifact rather than cohort size or biology."
        )
    elif surplus_survives:
        call = "data_property_not_sample_size"
        reading = (
            "N-matching leaves caudate's universe essentially unchanged while the surplus "
            "over both comparators survives, so the advantage is a property of the caudate "
            "data rather than of cohort size. Neither depth nor uniformity separates "
            "cleanly here -- inspect the paired contrasts before interpreting."
        )
    else:
        call = "no_material_surplus"
        reading = (
            "Neither sample size nor region explains a material universe difference once "
            "N is matched."
        )

    return [
        {
            "analysis_schema_version": ANALYSIS_SCHEMA_VERSION,
            "call": call,
            "reading": reading,
            "negligible_ratio_threshold": NEGLIGIBLE_RATIO,
            **verdict,
        }
    ]


def main() -> None:
    args = parse_args()
    indir = Path(args.indir)
    outdir = Path(args.outdir)

    passes = load_kind(indir, "coverage_pass")
    coverage = load_kind(indir, "donor_coverage")
    sd = load_kind(indir, "sd_cutoff")

    missing = sorted(set(passes.attrs["missing"]) | set(coverage.attrs["missing"]))
    if missing:
        msg = f"{len(missing)} region x chromosome outputs missing, e.g. {missing[:5]}"
        if args.require_complete:
            raise SystemExit(f"Incomplete: {msg}")
        print(f"WARNING: {msg}", file=sys.stderr)
        print("Results below are partial. Do not quote them.", file=sys.stderr)

    totals = genome_totals(passes)
    write_tsv(outdir / "universe_totals.tsv", totals.to_dict("records"))

    rows, verdict = summarize_rungs(totals, coverage)
    write_tsv(outdir / "universe_comparisons.tsv", rows)

    write_tsv(outdir / "sd_cutoff_summary.tsv", summarize_sd(sd))

    verdict["n_missing_region_chr"] = len(missing)
    verdict["autosomes_expected"] = len(AUTOSOMES)
    verdict_rows = build_verdict(verdict, totals)
    write_tsv(outdir / "universe_verdict.tsv", verdict_rows)

    # ---- console report -----------------------------------------------------
    print("\n=== Genome-wide CpG universe (autosomes) ===")
    full = totals[totals["design"] == "full_n"]
    for _, r in full.iterrows():
        print(
            f"  {r['region']:<12} N={DESIGN_N.get(r['region'], r['n_donors']):>3}  "
            f"{int(r['n_sites_pass']):>12,} pass  of {int(r['n_sites_pre']):>12,}  "
            f"({100 * r['n_sites_pass'] / r['n_sites_pre']:.1f}%)  [{int(r['n_chr'])} chr]"
        )

    matched = totals[(totals["design"] == "n_matched") & (totals["region"] == "caudate")]
    if not matched.empty:
        vals = sorted(matched["n_sites_pass"].tolist())
        print(
            f"  caudate      N=111  {statistics.fmean(vals):>12,.0f} pass  "
            f"(mean of {len(vals)} reps; range {vals[0]:,}-{vals[-1]:,})"
        )

    print("\n=== Per-donor coverage: depth vs uniformity ===")
    print(f"  {'region':<14}{'mean cov':>10}{'sd cov':>9}{'frac cov>=5':>13}{'sd frac':>9}{'min frac':>10}")
    for row in rows:
        if not row["comparison"].startswith("per_donor_coverage_"):
            continue
        print(
            f"  {row['comparison'].removeprefix('per_donor_coverage_'):<14}"
            f"{row['value_a']:>10}{row.get('sd_mean_cov', ''):>9}"
            f"{row['value_b']:>13}{row.get('sd_frac_ge_min', ''):>9}"
            f"{row.get('min_frac_ge_min', ''):>10}"
        )

    print("\n=== Comparison ladder ===")
    for row in rows:
        if row["comparison"].startswith("per_donor_coverage_"):
            continue
        ratio = f"ratio {row['ratio']}" if row["ratio"] != "" else ""
        print(f"  {row['comparison']:<44} {row['value_a']:>12}  {ratio}")

    print("\n=== Verdict ===")
    v = verdict_rows[0]
    print(f"  call: {v['call']}")
    print(f"  {v['reading']}")
    if missing:
        print(f"\n  PARTIAL: {len(missing)} of {len(REGIONS) * len(AUTOSOMES)} outputs missing.")


if __name__ == "__main__":
    main()
