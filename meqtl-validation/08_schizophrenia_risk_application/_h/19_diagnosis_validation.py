#!/usr/bin/env python3
"""Analysis 7: diagnosis association at prioritized SCZ-risk VMRs / linked features.

Orthogonal disease-layer check. Null case–control evidence does not
invalidate genetic meQTL evidence. Primary region: caudate AA.

Primary model (methylation / transcript feature):
  feature ~ primarydx + agedeath + sex + snpPC1–5

Primary + cell composition:
  + cellPC1–3 (MuSiC RNA)

Sensitivities (one at a time on top of primary):
  smoking, nicotine, cocaine, amphetamines, opioid (any of codeine/morphine/fentanyl),
  ethanol, antipsychotics, lifetime_antipsych, pmi, ph
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
MODULE = PROJECT / "meqtl-validation/08_schizophrenia_risk_application"
PHENOTYPE = PROJECT / "sample_summary/_m/phenotype_data.tsv"
SEED = 20260807

SENSITIVITY_VARS = [
    "smoking",
    "nicotine",
    "cocaine",
    "amphetamines",
    "opioid",
    "ethanol",
    "antipsychotics",
    "lifetime_antipsych",
    "pmi",
    "ph",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", default="caudate")
    p.add_argument("--population", default="AA", choices=["AA", "EA"])
    p.add_argument("--outdir", default="")
    p.add_argument("--fdr", type=float, default=0.05)
    return p.parse_args()


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    q_sorted = ranked * n / (np.arange(1, n + 1))
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)
    out = np.empty(n)
    out[order] = q_sorted
    q[ok] = out
    return q


def to_binary(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(float)
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    mapping = {
        "control": 0.0,
        "schizo": 1.0,
        "sczd": 1.0,
        "schizophrenia": 1.0,
        "m": 1.0,
        "f": 0.0,
        "male": 1.0,
        "female": 0.0,
        "true": 1.0,
        "false": 0.0,
        "yes": 1.0,
        "no": 0.0,
        "1": 1.0,
        "0": 0.0,
    }
    return series.astype(str).str.strip().str.lower().map(mapping)


def fit_ols(y: pd.Series, X: pd.DataFrame, predictor: str = "primarydx") -> dict:
    df = pd.concat([y.rename("y"), X], axis=1).dropna()
    if df.shape[0] < 20 or df["y"].nunique() < 2 or df[predictor].nunique() < 2:
        return {
            "n": int(df.shape[0]),
            "beta": np.nan,
            "se": np.nan,
            "t": np.nan,
            "pvalue": np.nan,
            "r2": np.nan,
        }
    model = sm.OLS(df["y"].astype(float), sm.add_constant(df[X.columns].astype(float), has_constant="add"))
    res = model.fit()
    return {
        "n": int(res.nobs),
        "beta": float(res.params.get(predictor, np.nan)),
        "se": float(res.bse.get(predictor, np.nan)),
        "t": float(res.tvalues.get(predictor, np.nan)),
        "pvalue": float(res.pvalues.get(predictor, np.nan)),
        "r2": float(res.rsquared),
    }


def load_prioritized_targets(module_m: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    pri = pd.read_csv(module_m / "prioritized/prioritized_loci.tsv", sep="\t")
    links = pd.read_csv(module_m / "level3/level3_vmr_feature_links.tsv", sep="\t")
    pri["locus_id"] = pri["locus_id"].astype(float)
    links["locus_id"] = links["locus_id"].astype(float)

    vmr_rows = []
    for _, row in pri.iterrows():
        vmrs = [str(int(float(row["best_vmr_id"])))]
        if pd.notna(row.get("tx_vmrs")) and str(row["tx_vmrs"]).strip():
            vmrs.extend([v.strip() for v in str(row["tx_vmrs"]).split(",") if v.strip()])
        for vid in sorted(set(vmrs)):
            vmr_rows.append(
                {
                    "locus_id": float(row["locus_id"]),
                    "index_snp": row["index_snp"],
                    "task_id": vid,
                    "is_best_vmr": vid == str(int(float(row["best_vmr_id"]))),
                }
            )
    vmr_df = pd.DataFrame(vmr_rows).drop_duplicates()

    sig = links[links["sig_tx_fdr"].astype(bool)].copy()
    sig = sig[sig["locus_id"].isin(pri["locus_id"])]
    # one best feature per locus×modality (lowest FDR)
    sig = sig.sort_values(["locus_id", "modality", "fdr"])
    feat_df = sig.groupby(["locus_id", "modality"], as_index=False).first()
    feat_df["task_id"] = feat_df["task_id"].astype(float).astype(int).astype(str)
    return vmr_df, feat_df


def load_covariate_frame(region: str, sample_ids: list[str]) -> pd.DataFrame:
    ph = pd.read_csv(PHENOTYPE, sep="\t")
    ph["brnum"] = ph["brnum"].astype(str)
    region_mask = ph["region"].astype(str).str.lower().str.replace(" ", "_")
    if region == "caudate":
        ph = ph[region_mask.str.contains("caudate")]
    elif region == "dlpfc":
        ph = ph[region_mask.str.contains("dlpfc|frontal")]
    else:
        ph = ph[region_mask.str.contains(region)]
    ph = ph[ph["brnum"].isin(sample_ids)].copy()
    ph = ph.drop_duplicates("brnum").set_index("brnum")

    out = pd.DataFrame(index=pd.Index(sample_ids, name="id"))
    out["primarydx"] = to_binary(ph.reindex(sample_ids)["primarydx"])
    out["agedeath"] = pd.to_numeric(ph.reindex(sample_ids)["agedeath"], errors="coerce")
    out["sex"] = to_binary(ph.reindex(sample_ids)["sex"])
    for pc in [f"snpPC{i}" for i in range(1, 6)]:
        out[pc] = pd.to_numeric(ph.reindex(sample_ids)[pc], errors="coerce")
    for col in [
        "smoking",
        "nicotine",
        "cocaine",
        "amphetamines",
        "ethanol",
        "antipsychotics",
        "lifetime_antipsych",
        "codeine",
        "morphine",
        "fentanyl",
        "pmi",
        "ph",
    ]:
        if col in ("pmi", "ph"):
            out[col] = pd.to_numeric(ph.reindex(sample_ids)[col], errors="coerce")
        else:
            out[col] = to_binary(ph.reindex(sample_ids)[col])
    opioid = out[["codeine", "morphine", "fentanyl"]].max(axis=1, skipna=True)
    # if all three missing, opioid missing
    all_miss = out[["codeine", "morphine", "fentanyl"]].isna().all(axis=1)
    out["opioid"] = opioid
    out.loc[all_miss, "opioid"] = np.nan

    cell_path = (
        PROJECT
        / "meqtl-validation/01_cpg_meqtl_mapping"
        / region
        / "_m/covariate_sensitivity/covariates_M5.txt"
    )
    if cell_path.exists():
        cell = pd.read_csv(cell_path, sep="\t")
        cell["id"] = cell["id"].astype(str)
        cell = cell.set_index("id")
        for c in ["cellPC1", "cellPC2", "cellPC3"]:
            if c in cell.columns:
                out[c] = pd.to_numeric(cell.reindex(sample_ids)[c], errors="coerce")
    return out


def load_vmr_mean_methylation(region: str, task_ids: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    prep = PROJECT / f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/prepared"
    bed_path = prep / "cpg_phenotypes.all_autosomes.bed.gz"
    if not bed_path.exists():
        raise SystemExit(f"Missing {bed_path}")

    wanted = set(str(x) for x in task_ids)
    # Discover chromosomes from elastic-net for requested task_ids
    en = pd.read_csv(
        PROJECT
        / "heritability/elastic_net_model/all_individuals"
        / region
        / "_m"
        / f"{region}_summary_elastic-net_AA.tsv",
        sep="\t",
    )
    en["task_id"] = en["task_id"].astype(str)
    en_sub = en[en["task_id"].isin(wanted)]
    chroms = sorted(en_sub["chrom"].astype(str).str.replace("^chr", "", regex=True).unique())

    cpg_to_vmr: dict[str, str] = {}
    for chrom in chroms:
        mpath = prep / f"cpg_vmr_map.chr{chrom}.tsv"
        if not mpath.exists():
            continue
        m = pd.read_csv(mpath, sep="\t", dtype={"vmr_id": str, "phenotype_id": str})
        m = m[m["vmr_id"].isin(wanted)]
        for _, r in m.iterrows():
            cpg_to_vmr[str(r["phenotype_id"])] = str(r["vmr_id"])

    if not cpg_to_vmr:
        raise SystemExit(f"No CpGs found for task_ids={sorted(wanted)}")

    with gzip.open(bed_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[4:]
        sums = {vid: np.zeros(len(samples), dtype=float) for vid in wanted}
        counts = {vid: np.zeros(len(samples), dtype=float) for vid in wanted}
        n_cpg = {vid: 0 for vid in wanted}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            pid = parts[3]
            if pid not in cpg_to_vmr:
                continue
            vid = cpg_to_vmr[pid]
            vals = np.array([np.nan if v in ("", "NA", "NaN") else float(v) for v in parts[4:]], dtype=float)
            ok = np.isfinite(vals)
            sums[vid][ok] += vals[ok]
            counts[vid][ok] += 1.0
            n_cpg[vid] += 1

    mat = {}
    meta_rows = []
    for vid in sorted(wanted):
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(counts[vid] > 0, sums[vid] / counts[vid], np.nan)
        mat[vid] = mean
        meta_rows.append(
            {
                "task_id": vid,
                "n_cpgs_in_mean": int(n_cpg[vid]),
                "mean_n_samples_per_cpg": float(np.nanmean(counts[vid])) if n_cpg[vid] else np.nan,
            }
        )
    y = pd.DataFrame(mat, index=samples)
    return y, pd.DataFrame(meta_rows)


def load_expression_features(feature_ids: list[str]) -> pd.DataFrame:
    bed = (
        PROJECT
        / "meqtl-validation/09_libd_eqtl_mapping/_m/caudate/genes/prepared/genes.expression.bed.gz"
    )
    if not bed.exists():
        raise SystemExit(f"Missing expression bed: {bed}")
    want_roots = {f.split(".")[0] for f in feature_ids}
    rows = {}
    with gzip.open(bed, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[4:]
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            pid = parts[3]
            root = pid.split(".")[0]
            if root in want_roots or pid in feature_ids:
                vals = [np.nan if v in ("", "NA", "NaN") else float(v) for v in parts[4:]]
                rows[pid] = vals
                # also store under requested id if version differs
                for f in feature_ids:
                    if f.split(".")[0] == root:
                        rows[f] = vals
    if not rows:
        return pd.DataFrame(index=samples)
    return pd.DataFrame(rows, index=samples)


def load_psi_features(feature_ids: list[str], region: str) -> pd.DataFrame:
    """Extract selected PSI features via R/SummarizedExperiment."""
    if not feature_ids:
        return pd.DataFrame()
    rse = PROJECT / f"inputs/counts/rse-psi.bsp3.{region}-n487.gencode-v47.RData"
    if not rse.exists():
        # try alternate naming
        cands = list((PROJECT / "inputs/counts").glob(f"rse-psi*{region}*"))
        if not cands:
            return pd.DataFrame()
        rse = cands[0]
    rbin = Path("/projects/p32505/opt/envs/rnaseq/bin/Rscript")
    if not rbin.exists():
        raise SystemExit(f"Rscript not found at {rbin}")

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        feat_file = tmp_path / "features.txt"
        out_file = tmp_path / "psi_matrix.tsv"
        feat_file.write_text("\n".join(feature_ids) + "\n")
        r_code = f"""
suppressPackageStartupMessages(library(SummarizedExperiment))
env <- new.env()
load('{rse}', envir=env)
se_names <- ls(envir=env)
se_names <- se_names[vapply(se_names, function(nm) inherits(env[[nm]], 'SummarizedExperiment'), logical(1))]
if (length(se_names) < 1) stop('No SummarizedExperiment in RData')
obj <- env[[se_names[1]]]
feats <- readLines('{feat_file}')
rn <- as.character(rownames(obj))
keep <- intersect(feats, rn)
if (length(keep) == 0) stop(paste('None of requested PSI features found. Example rownames:', paste(head(rn, 5), collapse=',')))
mat <- as.matrix(assay(obj)[keep, , drop=FALSE])
cd <- as.data.frame(colData(obj))
br <- if ('BrNum' %in% colnames(cd)) as.character(cd[['BrNum']]) else colnames(obj)
br[is.na(br) | br == ''] <- colnames(obj)[is.na(br) | br == '']
# Collapse technical replicates by BrNum
br_u <- unique(br)
out <- matrix(NA_real_, nrow=nrow(mat), ncol=length(br_u), dimnames=list(rownames(mat), br_u))
for (i in seq_along(br_u)) {{
  cols <- which(br == br_u[i])
  if (length(cols) == 1) out[, i] <- mat[, cols]
  else out[, i] <- rowMeans(mat[, cols, drop=FALSE], na.rm=TRUE)
}}
df <- data.frame(feature_id=rownames(out), out, check.names=FALSE)
write.table(df, '{out_file}', sep='\\t', quote=FALSE, row.names=FALSE)
"""
        script = tmp_path / "extract_psi.R"
        script.write_text(r_code)
        proc = subprocess.run([str(rbin), str(script)], check=False, capture_output=True, text=True)
        if proc.returncode != 0:
            sys.stderr.write(proc.stderr[-4000:] if proc.stderr else proc.stdout[-4000:])
            raise SystemExit(f"PSI extraction failed (exit {proc.returncode})")
        if not out_file.exists() or out_file.stat().st_size == 0:
            return pd.DataFrame()
        df = pd.read_csv(out_file, sep="\t")
        df = df.set_index("feature_id")
        return df.T


def model_matrix(cov: pd.DataFrame, model: str) -> pd.DataFrame:
    base = ["primarydx", "agedeath", "sex", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5"]
    cols = list(base)
    if model == "primary_cell":
        cols += ["cellPC1", "cellPC2", "cellPC3"]
    elif model.startswith("sens_"):
        var = model.replace("sens_", "", 1)
        cols.append(var)
    elif model == "primary_cell_plus_smoking":
        cols += ["cellPC1", "cellPC2", "cellPC3", "smoking"]
    keep = [c for c in cols if c in cov.columns]
    return cov[keep]


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir) if args.outdir else (MODULE / "_m/diagnosis")
    outdir.mkdir(parents=True, exist_ok=True)
    module_m = MODULE / "_m"

    vmr_df, feat_df = load_prioritized_targets(module_m)
    task_ids = sorted(vmr_df["task_id"].astype(str).unique())
    y_vmr, vmr_meta = load_vmr_mean_methylation(args.region, task_ids)
    sample_ids = list(y_vmr.index.astype(str))
    cov = load_covariate_frame(args.region, sample_ids)

    expr_ids = feat_df.loc[feat_df["modality"] == "expression", "feature_id"].astype(str).tolist()
    psi_ids = feat_df.loc[feat_df["modality"] == "psi", "feature_id"].astype(str).tolist()
    y_expr = load_expression_features(expr_ids) if expr_ids else pd.DataFrame()
    y_psi = load_psi_features(psi_ids, args.region) if psi_ids else pd.DataFrame()

    models = ["primary", "primary_cell"] + [f"sens_{v}" for v in SENSITIVITY_VARS] + [
        "primary_cell_plus_smoking"
    ]

    result_rows: list[dict] = []

    # VMR methylation tests
    for _, prow in vmr_df.iterrows():
        tid = str(prow["task_id"])
        if tid not in y_vmr.columns:
            continue
        y = y_vmr[tid]
        for model in models:
            X = model_matrix(cov, model)
            fit = fit_ols(y, X)
            result_rows.append(
                {
                    "feature_class": "vmr_mean_methylation",
                    "locus_id": prow["locus_id"],
                    "index_snp": prow["index_snp"],
                    "task_id": tid,
                    "feature_id": tid,
                    "gene_symbol": "",
                    "modality": "methylation",
                    "is_best_vmr": bool(prow["is_best_vmr"]),
                    "model": model,
                    **fit,
                }
            )

    # Transcript feature tests (align to available samples)
    for _, frow in feat_df.iterrows():
        fid = str(frow["feature_id"])
        modality = str(frow["modality"])
        src = y_expr if modality == "expression" else y_psi
        if src.empty or fid not in src.columns:
            # try versionless match for expression
            hit = None
            if modality == "expression" and not y_expr.empty:
                root = fid.split(".")[0]
                for c in y_expr.columns:
                    if c.split(".")[0] == root:
                        hit = c
                        break
            if hit is None:
                result_rows.append(
                    {
                        "feature_class": "linked_transcript",
                        "locus_id": frow["locus_id"],
                        "index_snp": frow["index_snp"],
                        "task_id": str(int(float(frow["task_id"]))),
                        "feature_id": fid,
                        "gene_symbol": frow.get("gene_symbol", ""),
                        "modality": modality,
                        "is_best_vmr": True,
                        "model": "primary",
                        "n": 0,
                        "beta": np.nan,
                        "se": np.nan,
                        "t": np.nan,
                        "pvalue": np.nan,
                        "r2": np.nan,
                        "note": "feature_matrix_missing",
                    }
                )
                continue
            fid_use = hit
        else:
            fid_use = fid
        y = src[fid_use]
        # covariates on overlapping BrNums
        common = sorted(set(y.index.astype(str)).intersection(cov.index.astype(str)))
        if not common:
            # expression bed may use BrNums not in meQTL set — rebuild cov for those IDs
            cov_tx = load_covariate_frame(args.region, list(y.index.astype(str)))
            common = sorted(set(y.index.astype(str)).intersection(cov_tx.index.astype(str)))
            cov_use = cov_tx
        else:
            cov_use = cov
        y_aligned = y.loc[common]
        cov_aligned = cov_use.loc[common]
        for model in models:
            X = model_matrix(cov_aligned, model)
            fit = fit_ols(y_aligned, X)
            result_rows.append(
                {
                    "feature_class": "linked_transcript",
                    "locus_id": frow["locus_id"],
                    "index_snp": frow["index_snp"],
                    "task_id": str(int(float(frow["task_id"]))),
                    "feature_id": fid,
                    "gene_symbol": frow.get("gene_symbol", ""),
                    "modality": modality,
                    "is_best_vmr": True,
                    "model": model,
                    **fit,
                }
            )

    res = pd.DataFrame(result_rows)
    # FDR within primary models only, separately by feature_class
    res["qvalue"] = np.nan
    for fclass in res["feature_class"].unique():
        mask = (res["feature_class"] == fclass) & (res["model"] == "primary")
        res.loc[mask, "qvalue"] = bh_fdr(res.loc[mask, "pvalue"].to_numpy())

    res_path = outdir / "diagnosis_association_results.tsv"
    res.to_csv(res_path, sep="\t", index=False)
    vmr_meta.to_csv(outdir / "vmr_mean_methylation_meta.tsv", sep="\t", index=False)
    vmr_df.to_csv(outdir / "tested_vmrs.tsv", sep="\t", index=False)
    feat_df.to_csv(outdir / "tested_transcript_features.tsv", sep="\t", index=False)

    # Snapshots
    primary_vmr = res[(res["feature_class"] == "vmr_mean_methylation") & (res["model"] == "primary")]
    primary_tx = res[(res["feature_class"] == "linked_transcript") & (res["model"] == "primary")]
    best_vmr_primary = primary_vmr[primary_vmr["is_best_vmr"].astype(bool)]

    def _n_sig(df: pd.DataFrame) -> int:
        if df.empty or "qvalue" not in df:
            return 0
        return int((pd.to_numeric(df["qvalue"], errors="coerce") < args.fdr).sum())

    # sensitivity: among primary-significant best VMRs, how many remain nominally p<0.05
    prim_sig = best_vmr_primary[pd.to_numeric(best_vmr_primary["qvalue"], errors="coerce") < args.fdr]
    sens_stability = []
    for var in SENSITIVITY_VARS + ["primary_cell", "primary_cell_plus_smoking"]:
        model = var if var.startswith("primary_") else f"sens_{var}"
        sub = res[
            (res["feature_class"] == "vmr_mean_methylation")
            & (res["is_best_vmr"].astype(bool))
            & (res["model"] == model)
            & (res["task_id"].isin(prim_sig["task_id"] if len(prim_sig) else []))
        ]
        if len(prim_sig) == 0:
            sens_stability.append(
                {
                    "model": model,
                    "n_primary_sig_best_vmrs": 0,
                    "n_remain_nominal_p0.05": 0,
                    "frac_remain": np.nan,
                }
            )
        else:
            n_remain = int((pd.to_numeric(sub["pvalue"], errors="coerce") < 0.05).sum())
            sens_stability.append(
                {
                    "model": model,
                    "n_primary_sig_best_vmrs": int(len(prim_sig)),
                    "n_remain_nominal_p0.05": n_remain,
                    "frac_remain": n_remain / len(prim_sig),
                }
            )
    sens_df = pd.DataFrame(sens_stability)
    sens_df.to_csv(outdir / "diagnosis_sensitivity_stability.tsv", sep="\t", index=False)

    n_best_vmr_sig = _n_sig(best_vmr_primary)
    n_any_vmr_sig = _n_sig(primary_vmr)
    n_tx_sig = _n_sig(primary_tx)
    any_dx = (n_any_vmr_sig + n_tx_sig) > 0
    # surviving smoking/antipsychotic among FDR-sig VMRs (best or secondary)
    prim_sig_any = primary_vmr[pd.to_numeric(primary_vmr["qvalue"], errors="coerce") < args.fdr]
    smoke_ok = ""
    ap_ok = ""
    if len(prim_sig_any) > 0:
        for label, model in [
            ("smoking", "sens_smoking"),
            ("antipsychotics", "sens_antipsychotics"),
            ("lifetime_antipsych", "sens_lifetime_antipsych"),
        ]:
            sub = res[
                (res["feature_class"] == "vmr_mean_methylation")
                & (res["model"] == model)
                & (res["task_id"].isin(prim_sig_any["task_id"]))
            ]
            n_remain = int((pd.to_numeric(sub["pvalue"], errors="coerce") < 0.05).sum())
            if label == "smoking":
                smoke_ok = bool(n_remain > 0)
            elif label == "antipsychotics":
                ap_ok = bool(n_remain > 0)
            else:
                ap_ok = bool(ap_ok) or bool(n_remain > 0)

    snapshot = {
        "region": args.region,
        "population": args.population,
        "n_prioritized_loci": int(vmr_df["locus_id"].nunique()),
        "n_vmrs_tested": int(vmr_df["task_id"].nunique()),
        "n_best_vmrs_tested": int(vmr_df.loc[vmr_df["is_best_vmr"], "task_id"].nunique()),
        "n_transcript_features_tested": int(feat_df.shape[0]),
        "n_best_vmr_primary_fdr_sig": n_best_vmr_sig,
        "n_any_vmr_primary_fdr_sig": n_any_vmr_sig,
        "n_transcript_primary_fdr_sig": n_tx_sig,
        "any_diagnosis_association_fdr": bool(any_dx),
        "primary_sig_survives_smoking_nominal": smoke_ok,
        "primary_sig_survives_antipsychotic_nominal": ap_ok,
        "interpretation": (
            "Diagnosis association detected among prioritized VMR features (FDR); supportive but not causal. "
            "Best-VMR-only FDR may be null while secondary TX-linked VMRs at the same loci are significant."
            if any_dx
            else (
                "No FDR-significant diagnosis association among prioritized VMR/transcript features; "
                "null does not invalidate genetic meQTL evidence."
            )
        ),
        "success_criterion": "optional_orthogonal_layer; null_allowed",
        "pass_optional_disease_layer": bool(any_dx),
    }
    write_tsv(outdir / "diagnosis_claim_snapshot.tsv", [snapshot])

    # Compact locus rollup (primary model)
    roll = []
    for locus_id, g in best_vmr_primary.groupby("locus_id"):
        g = g.sort_values("pvalue")
        top = g.iloc[0]
        tx = primary_tx[primary_tx["locus_id"] == locus_id]
        tx_top = tx.sort_values("pvalue").iloc[0] if len(tx) else None
        roll.append(
            {
                "locus_id": locus_id,
                "index_snp": top["index_snp"],
                "best_vmr_id": top["task_id"],
                "vmr_dx_beta": top["beta"],
                "vmr_dx_pvalue": top["pvalue"],
                "vmr_dx_qvalue": top["qvalue"],
                "vmr_dx_n": top["n"],
                "tx_feature_id": "" if tx_top is None else tx_top["feature_id"],
                "tx_modality": "" if tx_top is None else tx_top["modality"],
                "tx_gene_symbol": "" if tx_top is None else tx_top["gene_symbol"],
                "tx_dx_beta": "" if tx_top is None else tx_top["beta"],
                "tx_dx_pvalue": "" if tx_top is None else tx_top["pvalue"],
                "tx_dx_qvalue": "" if tx_top is None else tx_top["qvalue"],
                "tx_dx_n": "" if tx_top is None else tx_top["n"],
            }
        )
    write_tsv(outdir / "diagnosis_locus_rollup.tsv", roll)

    print(f"Wrote {res_path}")
    print(
        f"VMR FDR sig (best/any, primary): {n_best_vmr_sig}/{n_any_vmr_sig}; "
        f"transcript FDR sig: {n_tx_sig}"
    )


if __name__ == "__main__":
    main()
