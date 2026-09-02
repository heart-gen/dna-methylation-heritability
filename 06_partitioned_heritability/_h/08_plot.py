#!/usr/bin/env python3
"""06_partitioned_heritability -- figures for the partitioned-h2 result.

Usage:
    python 08_plot.py --run-id <id>

Two panels, both keyed on the tau z-score rather than the enrichment ratio:
enrichment on a continuous annotation is scale-dependent, so a bar chart of it
invites exactly the over-reading the module README warns against.

Brain traits and the prespecified non-brain controls are drawn together, because
the claim is comparative: "heritability concentrates in the annotation FOR BRAIN
TRAITS" is only supported if the controls do not behave the same way.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

CLASS_COLOR = {"brain": "#3B6FB6", "control": "#B0B0B0"}


def repo_root() -> Path:
    root = os.environ.get("V2_REPO_ROOT")
    if root:
        return Path(root)
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").is_dir():
            return parent
    raise SystemExit("Could not locate repository root")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "06_partitioned_heritability" / "_m" / "runs" / args.run_id
    metrics = run_dir / "results" / "sldsc-metrics.tsv"
    if not metrics.exists():
        raise SystemExit(f"No metrics table: {metrics} (run 07 first)")

    cfg = yaml.safe_load((root / "config" / "partitioned_heritability.yml").read_text())
    thr = float(cfg["fdr_threshold"])

    df = pd.read_csv(metrics, sep="\t").sort_values("tau_z")
    fig_dir = run_dir / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    ax = axes[0]
    colors = [CLASS_COLOR.get(c, "#888888") for c in df["trait_class"]]
    ypos = np.arange(len(df))
    ax.barh(ypos, df["tau_z"], color=colors)
    ax.set_yticks(ypos)
    ax.set_yticklabels(df["trait"])
    ax.axvline(0, color="black", lw=0.8)
    for z, label in ((1.96, "|z| = 1.96"), (-1.96, None)):
        ax.axvline(z, color="firebrick", ls="--", lw=0.8, label=label)
    ax.set_xlabel("Coefficient (tau) z-score")
    ax.set_title("Local SNP contribution annotation\nconditional on baselineLD")
    ax.legend(fontsize=8, loc="lower right")

    # Mark traits whose total h2 is not distinguishable from zero: their
    # enrichment is uninterpretable, and a reader should see that on the figure
    # rather than in a footnote.
    if "h2_interpretable" in df.columns:
        for i, ok in enumerate(df["h2_interpretable"]):
            if not bool(ok):
                ax.text(0, i, "  h2 n.s.", va="center", fontsize=7,
                        color="dimgray")

    ax = axes[1]
    q = df["tau_q"].astype(float).clip(lower=1e-300)
    ax.scatter(df["tau_z"], -np.log10(q), c=colors, s=45)
    for _, r in df.iterrows():
        ax.annotate(r["trait"], (r["tau_z"], -np.log10(max(float(r["tau_q"]), 1e-300))),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axhline(-np.log10(thr), color="firebrick", ls="--", lw=0.8,
               label=f"FDR {thr}")
    ax.set_xlabel("Coefficient (tau) z-score")
    ax.set_ylabel(r"$-\log_{10}$ FDR q")
    ax.set_title("Significance across the frozen trait family")
    ax.legend(fontsize=8)

    handles = [plt.Line2D([0], [0], marker="s", ls="", color=v,
                          label={"brain": "brain-relevant",
                                 "control": "non-brain control"}[k])
               for k, v in CLASS_COLOR.items()]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=9,
               frameon=False)

    run_label = df["run_id"].iloc[0] if "run_id" in df.columns else args.run_id
    fig.suptitle(f"Partitioned heritability -- {run_label}", fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    fig.savefig(fig_dir / "partitioned-h2.png", dpi=200)
    fig.savefig(fig_dir / "partitioned-h2.pdf")
    print(f"[06] wrote {fig_dir}/partitioned-h2.png")


if __name__ == "__main__":
    main()
