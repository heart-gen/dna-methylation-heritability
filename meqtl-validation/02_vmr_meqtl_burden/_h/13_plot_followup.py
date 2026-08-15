#!/usr/bin/env python3
"""Create the manuscript-facing calibrated-estimator follow-up figure."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REGIONS = ["caudate", "dlpfc", "hippocampus"]
REGION_LABEL = {"caudate": "Caudate", "dlpfc": "DLPFC", "hippocampus": "Hippocampus"}
REGION_COLOR = {"caudate": "#0072B2", "dlpfc": "#E69F00", "hippocampus": "#009E73"}
SENSITIVITY_COLOR = {
    "adjusted_prespecified": "#0072B2",
    "high_mappability": "#D55E00",
    "snp_proximity_excluded": "#009E73",
    "segdup_excluded": "#CC79A7",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-combined", required=True)
    parser.add_argument("--followup-combined", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def _tag(axis, label: str) -> None:
    axis.text(-0.15, 1.06, label, transform=axis.transAxes, fontsize=12, fontweight="bold")


def _forest(axis, frame: pd.DataFrame, y: np.ndarray, color: list[str], estimate: str, lower: str, upper: str) -> None:
    for position, (_, row) in zip(y, frame.iterrows()):
        point_color = color[int(position)] if int(position) < len(color) else "#444444"
        axis.errorbar(
            row[estimate], position,
            xerr=[[row[estimate] - row[lower]], [row[upper] - row[estimate]]],
            fmt="o", color=point_color, markersize=4, elinewidth=1, capsize=2,
        )


def main() -> None:
    args = parse_args()
    base = Path(args.base_combined)
    followup = Path(args.followup_combined)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)

    base_models = pd.read_csv(base / "burden_model_comparison.tsv", sep="\t")
    aligned = pd.read_csv(followup / "coordinate_aligned_models.tsv", sep="\t")
    boundary = pd.read_csv(followup / "boundary_excluded_models.tsv", sep="\t")
    annotation = pd.read_csv(followup / "annotation_models.tsv", sep="\t")
    external = pd.read_csv(followup / "external_validation_models.tsv", sep="\t")
    transcription = pd.read_csv(followup / "transcription_validation_models.tsv", sep="\t")

    plt.rcParams.update({
        "font.size": 8.5, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 7.5, "axes.linewidth": 0.8, "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.0), constrained_layout=True)

    # A: burden effect retained across coordinate and boundary definitions.
    ax = axes[0, 0]
    rows = []
    for region in REGIONS:
        primary = base_models[
            base_models["region"].eq(region) & base_models["test_family"].eq("calibrated_primary")
        ].iloc[0]
        coord = aligned[
            aligned["region"].eq(region) & aligned["test_family"].eq("calibrated_primary")
        ].iloc[0]
        bound = boundary[
            boundary["region"].eq(region)
            & boundary["vmr_definition"].eq("aa_exact_intersection")
            & boundary["test_family"].eq("calibrated_primary")
        ].iloc[0]
        for label, source in [
            ("Exact intersection", primary), ("Coordinate aligned", coord), ("Boundary excluded", bound)
        ]:
            rows.append({"region": region, "definition": label, **source.to_dict()})
    burden = pd.DataFrame(rows)
    definition_offset = {"Exact intersection": 0.18, "Coordinate aligned": 0.0, "Boundary excluded": -0.18}
    marker = {"Exact intersection": "o", "Coordinate aligned": "s", "Boundary excluded": "^"}
    for definition in definition_offset:
        d = burden[burden["definition"].eq(definition)]
        y = np.array([REGIONS[::-1].index(region) for region in d["region"]]) + definition_offset[definition]
        for index, (_, row) in enumerate(d.iterrows()):
            ax.errorbar(
                row["estimate_per_sd"], y[index],
                xerr=[[row["estimate_per_sd"] - row["ci_lower"]], [row["ci_upper"] - row["estimate_per_sd"]]],
                fmt=marker[definition], color=REGION_COLOR[row["region"]], markersize=4, capsize=2,
            )
    ax.axvline(0, color="#666666", linestyle="--", linewidth=0.7)
    ax.set_yticks(range(3), [REGION_LABEL[r] for r in REGIONS[::-1]])
    ax.set_xlabel("Change in meQTL-burden log odds per calibrated-estimate SD")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(
        handles=[Line2D([0], [0], marker=marker[d], color="#444444", linestyle="", label=d, markersize=5)
                 for d in definition_offset], frameon=False, loc="lower right",
    )
    _tag(ax, "A")

    # B: adjusted repeat/repressive effects.
    ax = axes[0, 1]
    feature_label = {"in_LINE_L1": "LINE/L1", "in_H3K9me3": "H3K9me3", "in_Quies": "Quiescent"}
    adj = annotation[
        annotation["estimator"].eq("calibrated")
        & annotation["sensitivity"].eq("adjusted_prespecified")
    ].copy()
    labels = []
    y = []
    colors = []
    for feature_index, feature in enumerate(feature_label):
        for region_index, region in enumerate(REGIONS):
            row = adj[adj["feature"].eq(feature) & adj["region"].eq(region)]
            if row.empty:
                continue
            labels.append(f"{feature_label[feature]} · {REGION_LABEL[region]}")
            y.append(len(labels) - 1)
            colors.append(REGION_COLOR[region])
    ordered = pd.concat([
        adj[adj["feature"].eq(feature) & adj["region"].eq(region)]
        for feature in feature_label for region in REGIONS
    ], ignore_index=True)
    _forest(ax, ordered, np.arange(len(ordered)), colors, "estimate_log_or_per_sd", "ci_lower", "ci_upper")
    ax.axvline(0, color="#666666", linestyle="--", linewidth=0.7)
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.invert_yaxis()
    ax.set_xlabel("Adjusted annotation log odds per calibrated-estimate SD")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    _tag(ax, "B")

    # C: technical restriction effects as a compact effect-size grid.
    ax = axes[1, 0]
    restrictions = list(SENSITIVITY_COLOR)
    technical = annotation[
        annotation["estimator"].eq("calibrated")
        & annotation["sensitivity"].isin(restrictions)
    ].copy()
    technical["row"] = technical["feature"].map(feature_label) + " · " + technical["region"].map(REGION_LABEL)
    row_order = [f"{feature_label[f]} · {REGION_LABEL[r]}" for f in feature_label for r in REGIONS]
    x = np.arange(len(restrictions))
    matrix = technical.pivot_table(index="row", columns="sensitivity", values="estimate_log_or_per_sd", aggfunc="first")
    matrix = matrix.reindex(index=row_order, columns=restrictions)
    vmax = float(np.nanmax(np.abs(matrix.to_numpy())))
    image = ax.imshow(matrix.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i, j]
            if pd.notna(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(value) > 0.65 * vmax else "#222222")
    ax.set_xticks(x, ["Adjusted", "High map.", "No SNP prox.", "No segdup"], rotation=20, ha="right")
    ax.set_yticks(np.arange(len(row_order)), row_order)
    fig.colorbar(image, ax=ax, fraction=0.045, pad=0.03, label="Log odds per SD")
    _tag(ax, "C")

    # D: independent external and existing transcriptional support.
    ax = axes[1, 1]
    ortho = []
    jaffe = external[
        external["estimator"].eq("calibrated") & external["sensitivity"].eq("primary")
    ]
    for _, row in jaffe.iterrows():
        ortho.append({"label": "Jaffe DLPFC meQTL", **row.to_dict()})
    tx = transcription[
        transcription["estimator"].eq("calibrated")
        & transcription["sensitivity"].eq("primary")
        & transcription["analysis"].isin(["transcription_expression", "transcription_psi"])
    ]
    for _, row in tx.iterrows():
        modality = "Expression" if row["analysis"].endswith("expression") else "Splicing"
        ortho.append({"label": f"{modality} · {REGION_LABEL[row['region']]}", **row.to_dict()})
    ortho_frame = pd.DataFrame(ortho)
    y = np.arange(len(ortho_frame))
    colors = [REGION_COLOR.get(row.get("region"), "#444444") for _, row in ortho_frame.iterrows()]
    _forest(ax, ortho_frame, y, colors, "estimate_log_or_per_sd", "ci_lower", "ci_upper")
    ax.axvline(0, color="#666666", linestyle="--", linewidth=0.7)
    x_left, x_right = ax.get_xlim()
    for row_index, row in ortho_frame.iterrows():
        if pd.isna(row.get("estimate_log_or_per_sd")):
            n_vmrs = pd.to_numeric(row.get("n_vmrs"), errors="coerce")
            n_label = f"; n={int(n_vmrs):,}" if pd.notna(n_vmrs) else ""
            ax.text(
                x_left + 0.03 * (x_right - x_left),
                row_index,
                f"Not estimable (no outcome variation{n_label})",
                ha="left",
                va="center",
                color="#555555",
                fontsize=7.5,
                fontstyle="italic",
            )
    ax.set_yticks(y, ortho_frame["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Association log odds per calibrated-estimate SD")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    _tag(ax, "D")

    for extension in ["pdf", "svg", "png"]:
        fig.savefig(output / f"calibrated-estimator-followup.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote follow-up figure under {output}")


if __name__ == "__main__":
    main()
