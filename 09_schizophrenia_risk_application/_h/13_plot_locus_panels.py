#!/usr/bin/env python3
"""Locus panels for the prioritized schizophrenia loci (09).

Usage:
    python _h/13_plot_locus_panels.py --run-id scz-AA-caudate-YYYYMMDD

One stacked regional panel per prioritized locus: GWAS association above, the
GTEx QTL track that colocalized best in the middle, the CpG meQTL track below,
with the PP4 of each ancestry-matched arm annotated. These are diagnostic
figures for reading a locus, not the manuscript figure -- Figure 5 is assembled
in 10_integrated_manuscript_outputs from this module's accepted tables.

The cross-ancestry meQTL arm is drawn because seeing the methylation signal is
the point of the panel, but its PP4 is labelled EXPLORATORY so a panel lifted
out of context cannot imply a colocalization claim it is not entitled to make.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def repo_root() -> Path:
    d = Path(__file__).resolve()
    while d != d.parent:
        if (d / ".git").is_dir():
            return d
        d = d.parent
    raise SystemExit("Could not locate repository root")


def neglog10(p: pd.Series) -> np.ndarray:
    # A p-value of exactly 0 underflows to inf and blows up the axis limits;
    # clip at the smallest positive double instead of dropping the variant.
    v = pd.to_numeric(p, errors="coerce").to_numpy(dtype=float)
    v = np.clip(v, np.nextafter(0, 1), 1.0)
    return -np.log10(v)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "09_schizophrenia_risk_application" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")
    fig_dir = run_dir / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    prio_f = run_dir / "results" / "prioritized-loci.tsv"
    prio = pd.read_csv(prio_f, sep="\t") if prio_f.is_file() else pd.DataFrame()
    if prio.empty or "prioritized" not in prio:
        (fig_dir / "no-panels.txt").write_text(
            "No prioritized locus; no locus panel was drawn.\n")
        print("No prioritized locus; nothing to plot")
        return
    prio = prio[prio["prioritized"].astype(bool)]

    abf_f = run_dir / "results" / "coloc-abf.tsv.gz"
    abf = pd.read_csv(abf_f, sep="\t") if abf_f.is_file() else pd.DataFrame()
    region_dir = run_dir / "coloc" / "regions"

    drawn = []
    for _, locus in prio.iterrows():
        locus_id = str(locus["locus_id"])
        rows = abf[(abf["locus_id"].astype(str) == locus_id)
                   & (abf["status"] == "OK")] if not abf.empty else pd.DataFrame()
        if rows.empty:
            print(f"locus {locus_id}: no evaluated coloc region; skipped")
            continue

        tracks = []
        matched = rows[rows["ld_ancestry_matched"].astype(bool)]
        if not matched.empty:
            tracks.append(matched.loc[matched["PP4"].idxmax()])
        unmatched = rows[~rows["ld_ancestry_matched"].astype(bool)]
        if not unmatched.empty:
            tracks.append(unmatched.loc[unmatched["PP4"].idxmax()])
        if not tracks:
            continue

        fig, axes = plt.subplots(len(tracks) + 1, 1, figsize=(9, 2.4 * (len(tracks) + 1)),
                                 sharex=True)
        axes = np.atleast_1d(axes)

        first = tracks[0]
        f = region_dir / (f"{first['chrom']}__{locus_id}__{first['arm']}__"
                          f"{first['qtl_context']}__"
                          f"{str(first['phenotype_id']).replace('/', '_')}.tsv.gz")
        if not f.is_file():
            plt.close(fig)
            continue
        gw = pd.read_csv(f, sep="\t")
        axes[0].scatter(gw["pos"] / 1e6, neglog10(gw["gwas_pvalue"]), s=6,
                        color="#4C4C4C", alpha=0.7, linewidths=0)
        axes[0].axhline(-np.log10(5e-8), color="#B22222", lw=0.8, ls="--")
        axes[0].set_ylabel("GWAS\n$-\\log_{10}P$")
        axes[0].set_title(
            f"Locus {locus_id} ({first['chrom']}) · index {locus.get('index_snp', 'NA')} · "
            f"PGC3 schizophrenia", fontsize=10, loc="left")

        for ax, tr in zip(axes[1:], tracks):
            tf = region_dir / (f"{tr['chrom']}__{locus_id}__{tr['arm']}__"
                               f"{tr['qtl_context']}__"
                               f"{str(tr['phenotype_id']).replace('/', '_')}.tsv.gz")
            if not tf.is_file():
                continue
            d = pd.read_csv(tf, sep="\t")
            matched_arm = bool(tr["ld_ancestry_matched"])
            colour = "#1F6FB4" if matched_arm else "#C77CFF"
            ax.scatter(d["pos"] / 1e6, neglog10(d["qtl_pvalue"]), s=6,
                       color=colour, alpha=0.7, linewidths=0)
            label = f"PP4 = {tr['PP4']:.3f}"
            if not matched_arm:
                label += "  (EXPLORATORY: cross-ancestry LD, not a coloc claim)"
            ax.set_ylabel(f"{tr['arm']}\n$-\\log_{{10}}P$")
            ax.set_title(f"{tr['qtl_context']} · {tr['phenotype_id']} · {label}",
                         fontsize=9, loc="left")

        axes[-1].set_xlabel(f"{first['chrom']} position (Mb, hg38)")
        fig.tight_layout()
        out = fig_dir / f"locus-{locus_id}-panel.png"
        fig.savefig(out, dpi=200)
        plt.close(fig)
        drawn.append({"locus_id": locus_id, "figure": out.name,
                      "n_tracks": len(tracks)})

    pd.DataFrame(drawn).to_csv(fig_dir / "locus-panels-index.tsv",
                               sep="\t", index=False)
    print(f"Drew {len(drawn)} locus panel(s) into {fig_dir}")


if __name__ == "__main__":
    main()
