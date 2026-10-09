#!/usr/bin/env python3
"""Targeted EA check of the illustrative loci, one region (09 stage 20).

Usage:
    python _h/20_ea_targeted_test.py --run-id scz-all_individuals.EA-crossregion-YYYYMMDD \
        --region dlpfc [--threads 4]

Tests, in the region's white-American (EA) donors, every risk-variant -> CpG
pair that stage 19 took from the AA primary, under the locked meQTL model M3a
re-estimated inside those donors (config/scz_ea_targeted.yml):

    CpG ~ dosage + agedeath + sex + primarydx + snpPC1-5 + methPC1-5

Steps, each written to work/<region>/ or results/<region>/:

  1. donors      the EA estimation cell's donor list (01b), by ID
  2. snpPC1-5    plink2 --pca 5 on those donors and the cell's own LD-pruned set
  3. identity    the pooled catalog's methylation for the AA donors must equal
                 the AA Module 05 run's tested methylation at every target CpG,
                 so the EA values are the phenotype the AA test used
  4. methPC1-5   the locked latent-factor recipe (seed, subsample, M0
                 residualization, SVD, sign rule), on the AA run's tested CpGs,
                 in the EA donors
  5. variants    REF/ALT must match the AA pfile per variant; EA MAF,
                 missingness and HWE under the locked genotype QC
  6. tests       tensorqtl's nominal statistic computed directly (residualize
                 dosage and phenotype on intercept + covariates; t from their
                 correlation; df = n - 2 - n_covariates); slope per ALT allele,
                 the AA convention
  7. cis scans   for the lead CpG of every locus whose lead pair reproduces,
                 the full +/- window scan the coloc stage reads
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

PLINK2 = "/projects/p32505/opt/bin/plink2"


def repo_root() -> Path:
    d = Path(__file__).resolve()
    while d != d.parent:
        if (d / ".git").is_dir():
            return d
        d = d.parent
    raise SystemExit("Could not locate repository root")


def sh(cmd: list[str]) -> None:
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f"command failed ({res.returncode}): {' '.join(cmd)}\n{res.stderr[-3000:]}")


# ------------------------------------------------------------------ methylation
class Matrix:
    """Column-subset reader for a Module 01 cpg_meth.phen (donors x every CpG).

    These files are up to ~2 million columns; only `cut` streams them safely
    (01_prepare_cpg_set.R records why fread cannot). Column names are bare
    positions; CpG IDs are chrN:pos.
    """

    def __init__(self, cat_dir: Path, chrom: str):
        self.f = cat_dir / "cpg" / f"chr_{chrom.replace('chr', '')}" / "cpg_meth.phen"
        if not self.f.is_file():
            raise SystemExit(f"missing {self.f}")
        with open(self.f) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
        if hdr[:2] != ["FID", "IID"]:
            raise SystemExit(f"{self.f}: unexpected header start {hdr[:2]}")
        self.chrom = chrom
        self.col = {f"{chrom}:{p}": i + 1 for i, p in enumerate(hdr) if i >= 2}

    def has(self, cpg_id: str) -> bool:
        return cpg_id in self.col

    def donors(self) -> list[str]:
        p = subprocess.run(["cut", "-f", "1", str(self.f)], capture_output=True,
                           text=True, check=True)
        return p.stdout.splitlines()[1:]

    def read(self, cpg_ids: list[str], donors: list[str]) -> pd.DataFrame:
        ids = [c for c in cpg_ids if c in self.col]
        if not ids:
            return pd.DataFrame(index=donors)
        fields = ",".join(str(x) for x in [1] + sorted(self.col[c] for c in ids))
        p = subprocess.run(["cut", "-f", fields, str(self.f)], capture_output=True,
                           text=True, check=True)
        lines = p.stdout.splitlines()
        hdr = lines[0].split("\t")
        df = pd.DataFrame([l.split("\t") for l in lines[1:]], columns=hdr)
        df = df.set_index("FID")
        df.columns = [f"{self.chrom}:{c}" for c in df.columns]
        df = df.apply(pd.to_numeric, errors="coerce")
        missing = [d for d in donors if d not in df.index]
        if missing:
            raise SystemExit(f"{self.f}: donors absent {missing[:5]}")
        return df.loc[donors, ids]


# ---------------------------------------------------------------- association
def residualize(M: np.ndarray, C: np.ndarray) -> np.ndarray:
    X = np.column_stack([np.ones(C.shape[0]), C])
    beta, *_ = np.linalg.lstsq(X, M, rcond=None)
    return M - X @ beta


def nominal(y: np.ndarray, g: np.ndarray, C: np.ndarray) -> dict:
    ok = np.isfinite(y)
    y, g, Cn = y[ok], g[ok], C[ok]
    g = np.where(np.isfinite(g), g, np.nanmean(g))
    yr = residualize(y[:, None], Cn)[:, 0]
    gr = residualize(g[:, None], Cn)[:, 0]
    n, k = len(y), Cn.shape[1]
    df = n - 2 - k
    if np.std(gr) == 0 or np.std(yr) == 0 or df < 1:
        return {"n": n, "df": df, "slope": np.nan, "slope_se": np.nan, "t": np.nan, "pval": np.nan}
    r = np.corrcoef(yr, gr)[0, 1]
    t = r * np.sqrt(df / (1 - r ** 2))
    slope = r * np.std(yr) / np.std(gr)
    return {"n": n, "df": df, "slope": slope, "slope_se": slope / t if t != 0 else np.nan,
            "t": t, "pval": 2 * stats.t.sf(abs(t), df)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--threads", default="4")
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "09_schizophrenia_gwas_loci" / "_m" / "runs" / args.run_id
    cfgd = run_dir / "code" / "config"
    cfg = yaml.safe_load((cfgd / "scz_ea_targeted.yml").read_text())
    covcfg = yaml.safe_load((cfgd / "covariates.yml").read_text())
    paths = yaml.safe_load((cfgd / "paths.yml").read_text())
    recipe = covcfg["primary_meqtl"]["latent_factor_recipe"]
    re_ = args.region
    work = run_dir / "work" / re_
    out = run_dir / "results" / re_
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    log = []

    targets = pd.read_csv(run_dir / "results" / "targets.tsv", sep="\t")
    targets = targets[targets["region"] == re_].copy()
    est = root / "01b_estimation_cells" / "_m" / "runs" / cfg["estimation_cells"][re_]
    cat = root / "01_vmr_catalog" / "_m" / "runs" / cfg["catalog_runs"][re_]
    cmb = root / "05_cpg_meqtl_burden" / "_m" / "runs" / cfg["meqtl_runs"][re_]

    # 1. donors ---------------------------------------------------------------
    keep = pd.read_csv(est / "vmr" / "donors_plink.txt", sep=r"\s+", header=None,
                       names=["FID", "IID"], dtype=str)
    if keep["FID"].duplicated().any():
        raise SystemExit("duplicate EA donor IDs")
    donors = keep["FID"].tolist()
    keep_f = work / "donors.keep"
    keep.to_csv(keep_f, sep="\t", header=False, index=False)
    log.append(("n_ea_donors", len(donors)))

    # 2. within-EA genotype PCs ----------------------------------------------
    pgen_prefix = str(root / paths["genotype"]["all_individuals"]["pgen"])[: -len(".pgen")]
    n_pc = int(cfg["ancestry_pcs"])
    sh([PLINK2, "--pfile", pgen_prefix, "--keep", str(keep_f), "--no-parents", "--no-sex",
        "--no-pheno", "--extract", str(est / "covs" / "pca" / "prune.prune.in"),
        "--pca", str(n_pc), "approx", "--threads", args.threads, "--out", str(work / "pca")])
    ev = pd.read_csv(work / "pca.eigenvec", sep=r"\s+")
    ev = ev.rename(columns={ev.columns[0]: "FID"}).set_index("FID")
    pcs = ev[[f"PC{i}" for i in range(1, n_pc + 1)]]
    pcs.columns = [f"snpPC{i}" for i in range(1, n_pc + 1)]
    pcs = pcs.loc[donors]

    # base covariates, from the estimation cell (01b subsets Module 01's files)
    covar = pd.read_csv(est / "covs" / "chr_1" / "TOPMed_LIBD.covar", sep=r"\s+", header=None,
                        names=["FID", "IID", "sex", "diagnosis"], dtype=str).set_index("FID")
    qcovar = pd.read_csv(est / "covs" / "chr_1" / "TOPMed_LIBD.qcovar", sep=r"\s+", header=None,
                         names=["FID", "IID", "agedeath"], dtype={0: str, 1: str}).set_index("FID")
    base = pd.concat([qcovar[["agedeath"]], covar[["sex", "diagnosis"]]], axis=1).loc[donors]
    base = pd.get_dummies(base, columns=["sex", "diagnosis"], drop_first=True, dtype=float)
    m0 = pd.concat([base, pcs], axis=1).astype(float)
    const = [c for c in m0.columns if m0[c].nunique() < 2]
    if const:
        raise SystemExit(f"constant covariate(s) in EA {re_}: {const}")

    # 3. identity: pooled-catalog methylation == AA tested methylation ---------
    id_rows = []
    for chrom, grp in targets.groupby("chrom"):
        mx = Matrix(cat, chrom)
        aa = pd.read_csv(cmb / "results" / "tested_meth" / f"{chrom}.tsv", sep="\t",
                         dtype={"FID": str}).set_index("FID")
        cps = [c for c in grp["cpg_id"].unique() if c in aa.columns and mx.has(c)]
        pooled_donors = set(mx.donors())
        shared = [d for d in aa.index if d in pooled_donors]
        if cps and shared:
            pooled = mx.read(cps, shared)
            diff = (pooled - aa.loc[shared, cps]).abs().to_numpy()
            id_rows.append((chrom, len(cps), len(shared), float(np.nanmax(diff))))
    ident = pd.DataFrame(id_rows, columns=["chrom", "n_cpgs", "n_aa_donors_shared", "max_abs_diff"])
    ident.to_csv(out / "phenotype-identity-check.tsv", sep="\t", index=False)
    if ident.empty or (ident["max_abs_diff"] > 1e-6).any():
        raise SystemExit(f"pooled-catalog methylation does not equal the AA tested methylation:\n{ident}")

    # 4. methPC1-5 in the EA donors -------------------------------------------
    mem = pd.read_csv(cmb / "results" / "tested-cpg-membership.tsv", sep="\t", dtype={"chr": str, "cpg_id": str})
    mem["chrom_n"] = mem["chr"].str.replace("chr", "", regex=False).astype(int)
    mem = mem.sort_values(["chrom_n", "cpg_pos", "cpg_id"]).reset_index(drop=True)
    mats = {c: Matrix(cat, c) for c in mem["chr"].unique()}
    ordered = [(c, i) for c, i in zip(mem["chr"], mem["cpg_id"]) if mats[c].has(i)]
    rng = np.random.default_rng(int(recipe["seed"]))
    n_sub = int(recipe["n_cpg_subsample"])
    sel = (np.sort(rng.choice(len(ordered), size=n_sub, replace=False))
           if len(ordered) > n_sub else np.arange(len(ordered)))
    chosen: dict[str, list[str]] = {}
    for i in sel:
        chosen.setdefault(ordered[i][0], []).append(ordered[i][1])
    Y = pd.concat([mats[c].read(ids, donors) for c, ids in chosen.items()], axis=1).to_numpy(float)
    sd = np.nanstd(Y, axis=0)
    Y = Y[:, np.isfinite(sd) & (sd > 0)]
    nmiss = int(np.isnan(Y).sum())
    if nmiss:
        cm = np.nanmean(Y, axis=0)
        idx = np.where(np.isnan(Y))
        Y[idx] = np.take(cm, idx[1])
    R = residualize(Y, m0.to_numpy())
    Z = (R - R.mean(0)) / np.clip(R.std(0), 1e-8, None)
    k = int(min(int(recipe["n_factors_estimated"]), Z.shape[0] - 1, Z.shape[1]))
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    for j in range(Vt.shape[0]):
        if Vt[j, np.argmax(np.abs(Vt[j]))] < 0:
            U[:, j] *= -1
            Vt[j] *= -1
    n_lf = len(covcfg["primary_meqtl"]["locked_latent_factors"])
    meth = pd.DataFrame((U[:, :k] * S[:k])[:, :n_lf], index=donors,
                        columns=[f"methPC{i}" for i in range(1, n_lf + 1)])
    cov = pd.concat([m0, meth], axis=1)
    cov.to_csv(out / "covariates.tsv", sep="\t")
    var_ratio = S ** 2 / np.sum(S ** 2)
    log += [("n_cpgs_latent_available", len(ordered)), ("n_cpgs_latent_used", Z.shape[1]),
            ("latent_seed", recipe["seed"]), ("n_missing_imputed", nmiss),
            ("var_explained_methPC1_5", float(var_ratio[:n_lf].sum())),
            ("covariates", ";".join(cov.columns)), ("n_covariates", cov.shape[1])]
    C = cov.to_numpy(float)

    # 5. variants: allele identity and EA QC ----------------------------------
    var_ids = sorted(targets["risk_variant_id"].unique())
    vfile = work / "target-variants.txt"
    vfile.write_text("\n".join(var_ids) + "\n")
    aa_prefix = str(root / paths["genotype"]["AA"]["pgen"])[: -len(".pgen")]

    def pvar_alleles(prefix: str) -> pd.DataFrame:
        p = subprocess.run(["grep", "-wFf", str(vfile), f"{prefix}.pvar"],
                           capture_output=True, text=True)
        rows = [l.split("\t")[:5] for l in p.stdout.splitlines() if not l.startswith("#")]
        return pd.DataFrame(rows, columns=["CHROM", "POS", "ID", "REF", "ALT"]).set_index("ID")

    pv_pool, pv_aa = pvar_alleles(pgen_prefix), pvar_alleles(aa_prefix)
    allele = pd.DataFrame(index=var_ids)
    allele["in_pooled_pfile"] = allele.index.isin(pv_pool.index)
    # The two pfiles share variant IDs but NOT always the REF/ALT assignment:
    # some variants are coded G/A in the AA pfile and A/G in the pooled one. The
    # AA slope is per AA-ALT allele, so a swapped variant is the same allele pair
    # counted the other way round and is harmonized by flipping the dosage. Only
    # a genuinely different allele set is excluded.
    def relation(v: str) -> str:
        if v not in pv_pool.index or v not in pv_aa.index:
            return "absent"
        pr, pa = pv_pool.loc[v, ["REF", "ALT"]], pv_aa.loc[v, ["REF", "ALT"]]
        if pr["REF"] == pa["REF"] and pr["ALT"] == pa["ALT"]:
            return "same"
        if pr["REF"] == pa["ALT"] and pr["ALT"] == pa["REF"]:
            return "swapped"
        return "different"
    allele["allele_relation_to_aa"] = [relation(v) for v in var_ids]
    allele["alleles_match_aa"] = allele["allele_relation_to_aa"].isin(["same", "swapped"])
    gq = cfg["genotype_qc"]
    qc_prefix = work / "variant-qc"
    sh([PLINK2, "--pfile", pgen_prefix, "--keep", str(keep_f), "--no-parents", "--no-sex",
        "--no-pheno", "--extract", str(vfile), "--freq", "--missing", "variant-only",
        "--hardy", "--threads", args.threads, "--out", str(qc_prefix)])
    fq = pd.read_csv(f"{qc_prefix}.afreq", sep=r"\s+").rename(columns={"#CHROM": "CHROM"}).set_index("ID")
    ms = pd.read_csv(f"{qc_prefix}.vmiss", sep=r"\s+").rename(columns={"#CHROM": "CHROM"}).set_index("ID")
    hw = pd.read_csv(f"{qc_prefix}.hardy", sep=r"\s+").rename(columns={"#CHROM": "CHROM"}).set_index("ID")
    allele["ea_alt_freq"] = fq["ALT_FREQS"].reindex(var_ids)
    allele["ea_maf"] = np.minimum(allele["ea_alt_freq"], 1 - allele["ea_alt_freq"])
    allele["ea_missing"] = ms["F_MISS"].reindex(var_ids)
    allele["ea_hwe_p"] = hw["P"].reindex(var_ids)
    allele["ea_qc_pass"] = (allele["alleles_match_aa"] & (allele["ea_maf"] >= gq["maf_min"]) &
                            (allele["ea_missing"] <= gq["missingness_max"]) &
                            (allele["ea_hwe_p"] >= gq["hwe_p_min"]))
    allele.index.name = "risk_variant_id"
    allele.to_csv(out / "variant-qc.tsv", sep="\t")

    def dosages(extract: Path, tag: str) -> pd.DataFrame:
        pre = work / f"dosage-{tag}"
        sh([PLINK2, "--pfile", pgen_prefix, "--keep", str(keep_f), "--no-parents", "--no-sex",
            "--no-pheno", "--extract", str(extract), "--export", "A", "--threads", args.threads,
            "--out", str(pre)])
        raw = pd.read_csv(f"{pre}.raw", sep=r"\s+", dtype={"FID": str})
        raw = raw.set_index("FID").drop(columns=["IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
        cols = {}
        for c in raw.columns:
            vid, counted = c.rsplit("_", 1)
            cols[c] = (vid, counted)
        out_df = pd.DataFrame({vid: raw[c].astype(float) for c, (vid, _) in cols.items()},
                              index=raw.index)
        out_df.attrs = {vid: counted for _, (vid, counted) in cols.items()}
        return out_df.loc[donors]

    # 6. tests ------------------------------------------------------------------
    D = dosages(vfile, "targets")
    pv_all = pv_pool
    for v in D.columns:
        counted = D.attrs.get(v)
        alt = pv_all.loc[v, "ALT"] if v in pv_all.index else None
        if counted != alt:
            if counted == (pv_all.loc[v, "REF"] if v in pv_all.index else None):
                D[v] = 2 - D[v]
            else:
                raise SystemExit(f"{v}: counted allele {counted} is neither REF nor ALT")
        # D now counts the pooled ALT; count the AA ALT instead where they differ.
        if v in allele.index and allele.loc[v, "allele_relation_to_aa"] == "swapped":
            D[v] = 2 - D[v]
    meth_t = pd.concat([Matrix(cat, c).read(sorted(g["cpg_id"].unique()), donors)
                        for c, g in targets.groupby("chrom")], axis=1)
    rows = []
    for _, r in targets.iterrows():
        base_row = r.to_dict()
        v, cpg = r["risk_variant_id"], r["cpg_id"]
        qc = allele.loc[v]
        status = "tested"
        if cpg not in meth_t.columns:
            status = "cpg_absent_in_pooled_matrix"
        elif not bool(qc["in_pooled_pfile"]) or v not in D.columns:
            status = "variant_absent_in_pooled_pfile"
        elif not bool(qc["alleles_match_aa"]):
            status = "alleles_differ_from_aa"
        elif not bool(qc["ea_qc_pass"]):
            status = "fails_ea_genotype_qc"
        res = {"n": np.nan, "df": np.nan, "slope": np.nan, "slope_se": np.nan, "t": np.nan, "pval": np.nan}
        if status == "tested":
            res = nominal(meth_t[cpg].to_numpy(float), D[v].to_numpy(float), C)
        base_row.update({f"ea_{k}": val for k, val in res.items()})
        base_row.update({"ea_status": status, "ea_maf": qc["ea_maf"], "ea_missing": qc["ea_missing"],
                         "ea_hwe_p": qc["ea_hwe_p"]})
        rows.append(base_row)
    res = pd.DataFrame(rows)
    res["sign_agrees"] = np.where(res["ea_status"] == "tested",
                                  np.sign(res["ea_slope"]) == np.sign(res["aa_slope"]), np.nan)
    res.to_csv(out / "pair-tests.tsv", sep="\t", index=False)

    # 7. cis scans for reproducing loci ---------------------------------------
    alpha = float(cfg["reading_rule"]["alpha"])
    win = int(cfg["coloc"]["window_bp"])
    lead = res[res["is_lead_pair"]]
    scans = []
    for _, r in lead.iterrows():
        if not (r["ea_status"] == "tested" and r["sign_agrees"] and r["ea_pval"] < alpha):
            continue
        cpg = r["cpg_id"]
        chrom, pos = cpg.split(":")
        pos = int(pos)
        pre = work / f"cis-{r['locus_id']}"
        sh([PLINK2, "--pfile", pgen_prefix, "--keep", str(keep_f), "--no-parents", "--no-sex",
            "--no-pheno", "--chr", chrom.replace("chr", ""), "--from-bp", str(max(1, pos - win)),
            "--to-bp", str(pos + win), "--snps-only", "just-acgt", "--maf", str(gq["maf_min"]),
            "--geno", str(gq["missingness_max"]), "--hwe", str(gq["hwe_p_min"]),
            "--write-snplist", "--threads", args.threads, "--out", str(pre)])
        Dc = dosages(Path(f"{pre}.snplist"), f"cis-{r['locus_id']}")
        ids = list(Dc.columns)
        y = meth_t[cpg].to_numpy(float) if cpg in meth_t.columns else None
        recs = []
        for v in ids:
            g = Dc[v].to_numpy(float)
            counted = Dc.attrs.get(v)
            vchrom, vpos, a1, a2 = v.split("_")[:4]
            nm = nominal(y, g, C)
            maf = np.nanmean(g) / 2
            recs.append({"variant_id": v, "chrom": f"chr{vchrom.replace('chr', '')}", "pos": int(vpos),
                         "counted_allele": counted, "id_allele_1": a1, "id_allele_2": a2,
                         "slope": nm["slope"], "slope_se": nm["slope_se"], "pval": nm["pval"],
                         "counted_freq": maf, "n": nm["n"]})
        sc = pd.DataFrame(recs)
        sc.insert(0, "cpg_id", cpg)
        sc.insert(0, "locus_id", r["locus_id"])
        sc.insert(0, "region", re_)
        sc.to_csv(out / f"cis-scan-locus{r['locus_id']}.tsv.gz", sep="\t", index=False)
        scans.append((r["locus_id"], cpg, len(sc), float(sc["pval"].min())))
    pd.DataFrame(scans, columns=["locus_id", "lead_cpg_id", "n_variants", "min_pval"]).to_csv(
        out / "cis-scan-summary.tsv", sep="\t", index=False)

    pd.DataFrame(log, columns=["field", "value"]).to_csv(out / "provenance.tsv", sep="\t", index=False)
    print(res.groupby("locus_id").agg(n=("ea_status", "size"),
                                      tested=("ea_status", lambda s: (s == "tested").sum()),
                                      lead_p=("ea_pval", "min")).to_string())


if __name__ == "__main__":
    main()
