#!/usr/bin/env python3
"""07_transcription_splicing_coupling -- figures.

Usage:
    python 05_plot.py --run-id <id>

Two panels: the coupling-test effect sizes per modality, and the coupled
fraction of VMRs across the local-genetic-control gradient. The second is the
descriptive form of the third test and is the one worth looking at when the
model is null but the gradient is not.
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

PRED_LABEL = {
    "any_meqtl_support": "any CpG meQTL support",
    "meqtl_proportion": "proportion meQTL-supported CpGs",
    "local_genetic_control": "local SNP contribution (z)",
}


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
    run_dir = root / "07_transcription_splicing_coupling" / "_m" / "runs" / args.run_id
    res_dir = run_dir / "results"
    tests_f = res_dir / "coupling-tests.tsv"
    if not tests_f.exists():
        raise SystemExit(f"No coupling tests: {tests_f}")
    tests = pd.read_csv(tests_f, sep="\t")

    fig_dir = res_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    ax = axes[0]
    t = tests[np.isfinite(tests["estimate"])].copy()
    if not t.empty:
        t["label"] = (t["modality"] + "\n" + t["predictor"].map(
            lambda p: PRED_LABEL.get(p, p)))
        t = t.sort_values("estimate")
        y = np.arange(len(t))
        lo = t["estimate"] - 1.96 * t["se"]
        hi = t["estimate"] + 1.96 * t["se"]
        sig = t["q"] < 0.05
        ax.errorbar(t["estimate"], y,
                    xerr=[t["estimate"] - lo, hi - t["estimate"]],
                    fmt="none", ecolor="#999999", lw=1)
        ax.scatter(t["estimate"], y,
                   c=["#C0392B" if s else "#3B6FB6" for s in sig], s=45,
                   zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(t["label"], fontsize=7)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlabel("log-odds of transcriptional coupling (95% CI)")
        ax.set_title("Coupling tests\n(red: FDR q < 0.05)")

    ax = axes[1]
    frame_f = res_dir / "coupling-model-frame.tsv"
    if frame_f.exists():
        fr = pd.read_csv(frame_f, sep="\t")
        for mod, g in fr.groupby("modality"):
            g = g[np.isfinite(g["local_snp_contribution_score_z"])]
            if g.empty:
                continue
            # Deciles of the continuous predictor: a descriptive view only.
            # The MODEL is fitted on the continuous score -- binning here is for
            # display and is never the tested form (AGENTS.md 3 bans grouping
            # the score as an analysis).
            try:
                q = pd.qcut(g["local_snp_contribution_score_z"], 10,
                            duplicates="drop")
            except ValueError:
                continue
            agg = g.groupby(q, observed=True).agg(
                frac=("coupled", "mean"),
                mid=("local_snp_contribution_score_z", "median"),
                n=("coupled", "size"))
            ax.plot(agg["mid"], agg["frac"], marker="o", ms=4, label=mod)
        ax.set_xlabel("local SNP contribution score (z), decile median")
        ax.set_ylabel("fraction of VMRs transcriptionally coupled")
        ax.set_title("Coupled fraction across the gradient\n(descriptive; model is continuous)")
        ax.legend(fontsize=8)

    fig.suptitle(f"Transcription / splicing coupling -- {args.run_id}", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(fig_dir / "coupling.png", dpi=200)
    fig.savefig(fig_dir / "coupling.pdf")
    print(f"[07] wrote {fig_dir}/coupling.png")


if __name__ == "__main__":
    main()
