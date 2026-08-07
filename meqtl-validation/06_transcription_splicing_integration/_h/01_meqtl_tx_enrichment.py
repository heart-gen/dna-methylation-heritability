#!/usr/bin/env python3
"""Phase 5: meQTL-supported / high-predictability VMRs vs expression and PSI links.

Reuses regulatory_context architecture_model_input tables (nearest-gene expression
and PSI window associations). Does not rerun transcriptome screens.

Primary tests:
  - any CpG meQTL support → transcriptional association enrichment
  - continuous predictability → transcriptional association
  - proportion meQTL-supported CpGs → transcriptional association
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
REGCTX = (
    PROJECT / "heritability/elastic_net_model/all_individuals/"
    "tissue_comparison/regulatory_context/_m"
)
PHASE2 = PROJECT / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m"
PHASE5 = PROJECT / "meqtl-validation" / "06_transcription_splicing_integration" / "_m"
SEED = 20260730

LINK_REL = {
    "expression": "{region}/AA/expression/nearest_gene_window_250kb/architecture_model_input.tsv",
    "psi": "{region}/AA/psi/window_250kb/architecture_model_input.tsv",
    "expression_abc": "{region}/AA/expression/abc/architecture_model_input.tsv",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True, choices=["caudate", "dlpfc", "hippocampus"])
    p.add_argument(
        "--modality",
        choices=["expression", "psi", "expression_abc"],
        default="expression",
    )
    p.add_argument("--burden-tsv", default="")
    p.add_argument("--link-tsv", default="")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=SEED)
    return p.parse_args()


def _z(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x * 0.0
    return (x - x.mean()) / sd


def arch_vmr_to_coord(vmr_id: str) -> str:
    """chr10_100932348_100932948 → 10:100932348-100932948"""
    parts = str(vmr_id).split("_")
    if len(parts) < 3:
        return str(vmr_id)
    chrom = parts[0].replace("chr", "")
    return f"{chrom}:{parts[1]}-{parts[2]}"


def load_links(region: str, modality: str, link_tsv: str) -> pd.DataFrame:
    path = Path(link_tsv) if link_tsv else REGCTX / LINK_REL[modality].format(region=region)
    if not path.exists():
        raise SystemExit(f"Missing link table: {path}")
    df = pd.read_csv(path, sep="\t")
    df["coord_id"] = df["vmr_id"].map(arch_vmr_to_coord)
    df["tx_associated"] = pd.to_numeric(df["any_sig_fdr_05"], errors="coerce").fillna(0).astype(int).eq(1)
    df["local_predictability_tx"] = pd.to_numeric(df["h2_unscaled"], errors="coerce")
    return df, path


def load_burden(region: str, burden_tsv: str) -> pd.DataFrame:
    path = Path(burden_tsv) if burden_tsv else PHASE2 / region / "vmr_meqtl_burden.tsv.gz"
    if not path.exists():
        raise SystemExit(f"Missing burden table: {path}")
    b = pd.read_csv(path, sep="\t")
    b["vmr_id"] = b["vmr_id"].astype(str)
    # Build coord keys from predictability summary for numeric task_ids
    pred = PROJECT / "heritability/elastic_net_model/all_individuals" / region / "_m" / f"{region}_summary_elastic-net_AA.tsv"
    if pred.exists():
        p = pd.read_csv(pred, sep="\t")
        p["task_id"] = p["task_id"].astype(str)
        p["coord_id"] = (
            p["chrom"].astype(str).str.replace("^chr", "", regex=True)
            + ":"
            + p["start"].astype(str)
            + "-"
            + p["end"].astype(str)
        )
        b = b.merge(
            p[["task_id", "coord_id"]].rename(columns={"task_id": "vmr_id"}),
            on="vmr_id",
            how="left",
        )
    # Also accept already-coordinate vmr_ids
    miss = b["coord_id"].isna() if "coord_id" in b.columns else pd.Series(True, index=b.index)
    if "coord_id" not in b.columns:
        b["coord_id"] = np.nan
    b.loc[miss, "coord_id"] = b.loc[miss, "vmr_id"].where(
        b.loc[miss, "vmr_id"].astype(str).str.contains(":"), np.nan
    )
    b["meqtl_supported"] = pd.to_numeric(b["n_cpgs_with_sig_meqtl"], errors="coerce").fillna(0).gt(0)
    return b


def fisher_enrichment(d: pd.DataFrame, xcol: str, ycol: str) -> dict:
    # Use 0/1 ints — pandas crosstab may coerce bool→uint8 and break False/True labels.
    x = d[xcol].fillna(False).astype(bool).astype(int)
    y = d[ycol].fillna(False).astype(bool).astype(int)
    tab = pd.crosstab(x, y).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
    or_, p = fisher_exact(tab.to_numpy(), alternative="greater")
    return {
        "odds_ratio": float(or_),
        "pvalue": float(p),
        "n00": int(tab.loc[0, 0]),
        "n01": int(tab.loc[0, 1]),
        "n10": int(tab.loc[1, 0]),
        "n11": int(tab.loc[1, 1]),
    }


def logistic_model(d: pd.DataFrame, ycol: str, xcols: list[str]) -> dict:
    use = d.dropna(subset=[ycol] + xcols).copy()
    if use.empty or use[ycol].nunique() < 2:
        return {"coef": np.nan, "pvalue": np.nan, "n": len(use), "error": "insufficient"}
    y = use[ycol].astype(float)
    X = use[xcols].apply(lambda s: pd.to_numeric(s, errors="coerce"))
    # Binary predictors: keep 0/1; continuous predictors: z-score
    for c in xcols:
        if X[c].dropna().nunique() <= 2:
            X[c] = X[c].fillna(0).astype(float)
        else:
            X[c] = _z(X[c])
    X = sm.add_constant(X, has_constant="add")
    try:
        res = sm.GLM(y, X, family=sm.families.Binomial()).fit()
        primary = xcols[0]
        return {
            "coef": float(res.params[primary]),
            "or": float(np.exp(res.params[primary])),
            "pvalue": float(res.pvalues[primary]),
            "n": int(len(use)),
            "covariates": ",".join(xcols),
            "error": "",
        }
    except Exception as exc:  # noqa: BLE001
        return {"coef": np.nan, "pvalue": np.nan, "n": len(use), "error": str(exc)}


def main() -> None:
    args = parse_args()
    region = args.region
    links, link_path = load_links(region, args.modality, args.link_tsv)
    burden = load_burden(region, args.burden_tsv)

    # Join on coordinates
    bkeep = [
        c for c in [
            "coord_id", "vmr_id", "meqtl_supported", "n_tested_cpgs",
            "n_cpgs_with_sig_meqtl", "proportion_cpgs_with_sig_meqtl",
            "local_predictability", "average_cpg_coverage", "mean_cpg_variance",
            "length", "umap_k24_mean",
        ] if c in burden.columns
    ]
    d = links.merge(burden[bkeep], on="coord_id", how="inner", suffixes=("_tx", "_meqtl"))
    # Prefer Phase 2 predictability when present
    if "local_predictability" in d.columns:
        d["predictability"] = pd.to_numeric(d["local_predictability"], errors="coerce")
    else:
        d["predictability"] = d["local_predictability_tx"]
    d["meqtl_supported"] = d["meqtl_supported"].fillna(False).astype(bool)
    d["tx_associated"] = d["tx_associated"].astype(bool)
    d["both_tx"] = d["tx_associated"]  # placeholder; filled in summary for expression∩psi later

    if d.empty:
        raise SystemExit(
            f"No overlapping VMRs between burden and {args.modality} links for {region}"
        )

    rows = []
    # 1) Fisher: meQTL support × tx association
    fish = fisher_enrichment(d, "meqtl_supported", "tx_associated")
    rows.append({
        "region": region,
        "modality": args.modality,
        "model": "fisher_meqtl_support_x_tx",
        "n_vmrs": len(d),
        "n_meqtl_supported": int(d["meqtl_supported"].sum()),
        "n_tx_associated": int(d["tx_associated"].sum()),
        "n_both": fish["n11"],
        "estimate": fish["odds_ratio"],
        "pvalue": fish["pvalue"],
        "covariates": "",
        "link_path": str(link_path),
    })

    # 2) Logistic: tx ~ meQTL support + technical
    adj = []
    for c in ["n_features_tested", "vmr_length", "min_distance", "methylation_variance", "num_snps"]:
        if c in d.columns and d[c].notna().sum() >= 50:
            adj.append(c)
    # map length if needed
    if "vmr_length" not in d.columns and "length" in d.columns:
        d["vmr_length"] = d["length"]
        if "vmr_length" not in adj:
            adj.append("vmr_length")
    d["y"] = d["tx_associated"].astype(int)
    d["meqtl_supported_int"] = d["meqtl_supported"].astype(int)
    log1 = logistic_model(d, "y", ["meqtl_supported_int"] + adj)
    rows.append({
        "region": region,
        "modality": args.modality,
        "model": "logistic_tx_~_meqtl_support_adjusted",
        "n_vmrs": log1.get("n", len(d)),
        "n_meqtl_supported": int(d["meqtl_supported"].sum()),
        "n_tx_associated": int(d["tx_associated"].sum()),
        "n_both": int((d["meqtl_supported"] & d["tx_associated"]).sum()),
        "estimate": log1.get("or", log1.get("coef")),
        "pvalue": log1.get("pvalue"),
        "covariates": log1.get("covariates", ""),
        "link_path": str(link_path),
        "error": log1.get("error", ""),
    })

    # 3) Logistic: tx ~ continuous predictability
    log2 = logistic_model(d, "y", ["predictability"] + adj)
    rows.append({
        "region": region,
        "modality": args.modality,
        "model": "logistic_tx_~_predictability_adjusted",
        "n_vmrs": log2.get("n", np.nan),
        "n_meqtl_supported": int(d["meqtl_supported"].sum()),
        "n_tx_associated": int(d["tx_associated"].sum()),
        "n_both": int((d["meqtl_supported"] & d["tx_associated"]).sum()),
        "estimate": log2.get("or", log2.get("coef")),
        "pvalue": log2.get("pvalue"),
        "covariates": log2.get("covariates", ""),
        "link_path": str(link_path),
        "error": log2.get("error", ""),
    })

    # 4) Logistic: tx ~ proportion meQTL CpGs
    if "proportion_cpgs_with_sig_meqtl" in d.columns:
        log3 = logistic_model(d, "y", ["proportion_cpgs_with_sig_meqtl"] + adj)
        rows.append({
            "region": region,
            "modality": args.modality,
            "model": "logistic_tx_~_meqtl_proportion_adjusted",
            "n_vmrs": log3.get("n", np.nan),
            "n_meqtl_supported": int(d["meqtl_supported"].sum()),
            "n_tx_associated": int(d["tx_associated"].sum()),
            "n_both": int((d["meqtl_supported"] & d["tx_associated"]).sum()),
            "estimate": log3.get("or", log3.get("coef")),
            "pvalue": log3.get("pvalue"),
            "covariates": log3.get("covariates", ""),
            "link_path": str(link_path),
            "error": log3.get("error", ""),
        })

    # 5) Mann-Whitney: predictability in tx-associated vs not
    a = d.loc[d["tx_associated"], "predictability"].dropna()
    b = d.loc[~d["tx_associated"], "predictability"].dropna()
    if len(a) and len(b):
        stat, p = mannwhitneyu(a, b, alternative="greater")
        rows.append({
            "region": region,
            "modality": args.modality,
            "model": "mannwhitney_predictability_tx_greater",
            "n_vmrs": len(d),
            "n_meqtl_supported": int(d["meqtl_supported"].sum()),
            "n_tx_associated": int(d["tx_associated"].sum()),
            "n_both": int((d["meqtl_supported"] & d["tx_associated"]).sum()),
            "estimate": float(a.median() - b.median()),
            "pvalue": float(p),
            "covariates": "",
            "link_path": str(link_path),
            "median_pred_tx": float(a.median()),
            "median_pred_nontx": float(b.median()),
        })

    outdir = PHASE5 / region
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"meqtl_x_{args.modality}_enrichment.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)

    # VMR-level joined table for downstream locus picks
    keep = [
        c for c in [
            "coord_id", "vmr_id_tx", "vmr_id_meqtl", "vmr_id", "seqnames", "start", "end",
            "tx_associated", "meqtl_supported", "predictability",
            "proportion_cpgs_with_sig_meqtl", "n_sig_fdr_05", "n_features_tested",
            "min_distance", "max_abs_beta", "chromatin_state", "h2_category",
        ] if c in d.columns or c in {"vmr_id_tx", "vmr_id_meqtl"}
    ]
    # normalize vmr id columns after merge suffixes
    if "vmr_id_tx" not in d.columns and "vmr_id" in links.columns:
        pass
    d_out = d.copy()
    if "vmr_id_x" in d_out.columns:
        d_out = d_out.rename(columns={"vmr_id_x": "vmr_id_tx", "vmr_id_y": "vmr_id_meqtl"})
    d_out.to_csv(outdir / f"vmr_meqtl_{args.modality}_joined.tsv.gz", sep="\t", index=False, compression="gzip")
    print(
        f"{region} {args.modality}: n={len(d)} meQTL+={d['meqtl_supported'].sum()} "
        f"tx+={d['tx_associated'].sum()} both={(d['meqtl_supported'] & d['tx_associated']).sum()} "
        f"Fisher OR={fish['odds_ratio']:.3f} p={fish['pvalue']:.2e} → {out}"
    )


if __name__ == "__main__":
    main()
