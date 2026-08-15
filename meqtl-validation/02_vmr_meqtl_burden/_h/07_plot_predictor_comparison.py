#!/usr/bin/env python3
"""Create the manuscript-facing calibrated-versus-legacy comparison figure."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REGIONS = ["caudate", "dlpfc", "hippocampus"]
REGION_LABELS = {"caudate": "Caudate", "dlpfc": "DLPFC", "hippocampus": "Hippocampus"}
REGION_COLORS = {"caudate": "#0072B2", "dlpfc": "#E69F00", "hippocampus": "#009E73"}
PREDICTOR_COLORS = {"calibrated": "#0072B2", "legacy": "#D55E00"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-root", required=True)
    parser.add_argument("--combined-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seed", type=int, default=20260808)
    return parser.parse_args()


def _panel_tag(axis, label: str) -> None:
    axis.text(-0.14, 1.06, label, transform=axis.transAxes, fontsize=12, fontweight="bold")


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def main() -> None:
    args = parse_args()
    region_root = Path(args.region_root)
    combined = Path(args.combined_dir)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    bridges = []
    for region in REGIONS:
        frame = pd.read_csv(region_root / region / "predictor_bridge.tsv.gz", sep="\t")
        frame = frame[_as_bool(frame["common_complete_case"])].copy()
        if len(frame) > 2500:
            frame = frame.iloc[np.sort(rng.choice(len(frame), 2500, replace=False))]
        frame["region"] = region
        bridges.append(frame)
    bridge = pd.concat(bridges, ignore_index=True)
    models = pd.read_csv(combined / "burden_model_comparison.tsv", sep="\t")
    matched = pd.read_csv(combined / "matched_comparison.tsv", sep="\t")
    transition = pd.read_csv(combined / "quintile_transition.tsv", sep="\t")

    plt.rcParams.update({
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.6), constrained_layout=True)

    ax = axes[0, 0]
    for region in REGIONS:
        data = bridge[bridge["region"].eq(region)]
        ax.scatter(
            data["legacy_h2_unscaled"], data["calibrated_h2"],
            s=5, alpha=0.10, linewidths=0, color=REGION_COLORS[region], rasterized=True,
        )
    ax.axhline(0, color="#666666", linewidth=0.6, linestyle="--")
    ax.set_xlabel("Legacy aggregate SNP contribution")
    ax.set_ylabel("Simulation-calibrated local SNP variance")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        handles=[Line2D([0], [0], marker="o", linestyle="", color=REGION_COLORS[r],
                        label=REGION_LABELS[r], markersize=5) for r in REGIONS],
        frameon=False, loc="upper left",
    )
    _panel_tag(ax, "A")

    ax = axes[0, 1]
    matrix = transition.groupby(
        ["legacy_quintile", "calibrated_quintile"], observed=False
    )["n_vmrs"].sum().unstack(fill_value=0)
    matrix = matrix.reindex(index=[f"Q{i}" for i in range(1, 6)], columns=[f"Q{i}" for i in range(1, 6)])
    row_fraction = matrix.div(matrix.sum(axis=1), axis=0)
    image = ax.imshow(row_fraction.to_numpy(), cmap="Blues", vmin=0, vmax=max(0.2, row_fraction.max().max()))
    for i in range(5):
        for j in range(5):
            value = row_fraction.iloc[i, j]
            ax.text(j, i, f"{100 * value:.0f}%", ha="center", va="center",
                    color="white" if value > 0.35 else "#222222", fontsize=8)
    ax.set_xticks(range(5), row_fraction.columns)
    ax.set_yticks(range(5), row_fraction.index)
    ax.set_xlabel("Calibrated-estimate quintile")
    ax.set_ylabel("Legacy-estimate quintile")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03, label="Row fraction")
    _panel_tag(ax, "B")

    ax = axes[1, 0]
    forest = models[models["test_family"].isin(["calibrated_primary", "legacy_comparator"])].copy()
    y_positions = {region: i for i, region in enumerate(REGIONS[::-1])}
    offsets = {"calibrated_primary": 0.12, "legacy_comparator": -0.12}
    labels = {"calibrated_primary": "Calibrated", "legacy_comparator": "Legacy"}
    colors = {"calibrated_primary": PREDICTOR_COLORS["calibrated"],
              "legacy_comparator": PREDICTOR_COLORS["legacy"]}
    for family in ["calibrated_primary", "legacy_comparator"]:
        data = forest[forest["test_family"].eq(family)]
        for _, row in data.iterrows():
            y = y_positions[row["region"]] + offsets[family]
            ax.errorbar(
                row["estimate_per_sd"], y,
                xerr=[[row["estimate_per_sd"] - row["ci_lower"]],
                      [row["ci_upper"] - row["estimate_per_sd"]]],
                fmt="o", color=colors[family], capsize=2, markersize=4, linewidth=1,
            )
    ax.axvline(0, color="#666666", linewidth=0.7, linestyle="--")
    ax.set_yticks([y_positions[r] for r in REGIONS[::-1]], [REGION_LABELS[r] for r in REGIONS[::-1]])
    ax.set_xlabel("Change in meQTL-burden log odds per predictor SD")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(
        handles=[Line2D([0], [0], marker="o", linestyle="", color=colors[f],
                        label=labels[f], markersize=5)
                 for f in ["calibrated_primary", "legacy_comparator"]],
        frameon=False, loc="lower right",
    )
    _panel_tag(ax, "C")

    ax = axes[1, 1]
    match_order = ["calibrated_positive_vs_nonpositive", "legacy_high_vs_low"]
    match_labels = {match_order[0]: "Calibrated positive", match_order[1]: "Legacy high"}
    match_colors = {match_order[0]: PREDICTOR_COLORS["calibrated"],
                    match_order[1]: PREDICTOR_COLORS["legacy"]}
    x = np.arange(len(REGIONS))
    width = 0.34
    for index, analysis in enumerate(match_order):
        data = matched[matched["analysis"].eq(analysis)].set_index("region").reindex(REGIONS)
        ax.bar(
            x + (index - 0.5) * width,
            data["mean_difference"], width=width,
            color=match_colors[analysis], label=match_labels[analysis], alpha=0.9,
        )
    ax.axhline(0, color="#666666", linewidth=0.7)
    ax.set_xticks(x, [REGION_LABELS[r] for r in REGIONS], rotation=15, ha="right")
    ax.set_ylabel("Matched difference in meQTL burden")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper right")
    _panel_tag(ax, "D")

    for extension in ["pdf", "svg", "png"]:
        fig.savefig(output / f"predictor-comparison.{extension}", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote comparison figure under {output}")


if __name__ == "__main__":
    main()
