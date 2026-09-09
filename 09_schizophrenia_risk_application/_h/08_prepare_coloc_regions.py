#!/usr/bin/env python3
"""Assemble harmonised per-locus regions for colocalization (09).

Usage (inside an array task):
    python _h/08_prepare_coloc_regions.py --run-id scz-AA-caudate-YYYYMMDD --chrom 22

For every schizophrenia locus on one chromosome this writes one long-format
table per (arm, phenotype) holding the GWAS and QTL statistics on a COMMON set
of variants:

    arm_a  PGC3 EUR GWAS x GTEx v11 brain eQTL / sQTL   (EUR x EUR)
    arm_b  PGC3 EUR GWAS x this cohort's CpG meQTL      (EUR x AA)

Harmonisation is the whole job, and it is where colocalization silently goes
wrong. Three rules, applied here so no later stage can skip them:

  * variants are matched on hg38 position AND the unordered allele pair, never
    on position alone -- position-only matching joins multi-allelic sites to the
    wrong allele and produces confident nonsense;
  * when a study reports the other allele as its effect allele, the effect SIGN
    is flipped rather than the variant dropped;
  * strand-ambiguous variants (A/T and C/G) are DROPPED, because the two studies
    were not necessarily strand-aligned and an unflippable ambiguity cannot be
    resolved from summary statistics.

`n_variants_shared` is recorded per region. A thin overlap does not make a
region "not colocalised"; it makes it unevaluable, and stage 09 records it that
way rather than reporting a posterior computed from a handful of variants.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.compute as pc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_config import load_run_config  # noqa: E402

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}


def repo_root() -> Path:
    d = Path(__file__).resolve()
    while d != d.parent:
        if (d / ".git").is_dir():
            return d
        d = d.parent
    raise SystemExit("Could not locate repository root")


def read_manifest(run_dir: Path) -> dict:
    m = pd.read_csv(run_dir / "manifest.tsv", sep="\t", dtype=str)
    return dict(zip(m["field"], m["value"]))


def split_variant(variant_id: pd.Series) -> pd.DataFrame:
    """chrom_pos_ref_alt[_...] -> components. Non-conforming rows become NaN."""
    parts = variant_id.astype(str).str.split("_", expand=True)
    if parts.shape[1] < 4:
        raise SystemExit("Variant IDs do not have chrom_pos_ref_alt structure")
    return pd.DataFrame({
        "chrom": parts[0],
        "pos": pd.to_numeric(parts[1], errors="coerce"),
        "a1": parts[2].str.upper(),
        "a2": parts[3].str.upper(),
    })


def add_keys(df: pd.DataFrame) -> pd.DataFrame:
    """Unordered-allele join key, plus the ambiguity flag."""
    df = df.copy()
    pair = [frozenset((x, y)) for x, y in zip(df["a1"], df["a2"])]
    df["strand_ambiguous"] = [p in AMBIGUOUS for p in pair]
    df["match_key"] = [
        f"{c}_{int(p)}_{'_'.join(sorted((x, y)))}"
        for c, p, x, y in zip(df["chrom"], df["pos"], df["a1"], df["a2"])
    ]
    return df


def harmonise(gwas: pd.DataFrame, qtl: pd.DataFrame) -> pd.DataFrame:
    """Inner-join on match_key and orient the QTL effect to the GWAS allele.

    Only biallelic single-nucleotide alleles are keyed; indels whose allele
    strings differ between resources would not match anyway, and silently
    keeping them would mean an inner join that looks clean but drops them
    without a count.
    """
    g = gwas[~gwas["strand_ambiguous"]]
    q = qtl[~qtl["strand_ambiguous"]]
    # A duplicated key means the same position/allele pair appears twice in one
    # resource; there is no basis for choosing between them, so drop both.
    g = g[~g["match_key"].duplicated(keep=False)]
    q = q[~q["match_key"].duplicated(keep=False)]

    m = g.merge(q, on="match_key", suffixes=("_gwas", "_qtl"))
    if m.empty:
        return m
    same = m["a1_gwas"].values == m["a1_qtl"].values
    flipped = m["a1_gwas"].values == np.array(
        [COMPLEMENT.get(a, "N") for a in m["a1_qtl"].values])
    aligned = same | flipped
    m["qtl_beta_aligned"] = np.where(aligned, m["beta_qtl"], -m["beta_qtl"])
    m["allele_orientation"] = np.where(aligned, "same", "swapped")
    return m


# GTEx v11 does not name the two all-pairs resources alike: eQTL files are
# {tissue}.v11.allpairs.{chrom}.parquet, sQTL files are
# {tissue}.v11.cis_sqtl.allpairs.{chrom}.parquet. Both are keyed here so a
# pattern mismatch cannot silently drop an entire arm.
ALLPAIRS_PATTERN = {
    "gtex_eqtl": "{tissue}.v11.allpairs.{chrom}.parquet",
    "gtex_sqtl": "{tissue}.v11.cis_sqtl.allpairs.{chrom}.parquet",
}


def in_any_window(pos: pd.Series, windows: list[tuple[int, int]]) -> pd.Series:
    """Row mask for positions inside any locus window.

    Filtering on the min-max SPAN instead would keep nearly the whole
    chromosome whenever a chromosome's loci are far apart -- on chr22 that is
    ~50 Mb rather than the ~7 Mb the three windows actually cover, and it is
    what made a single sQTL tissue too large to hold.
    """
    mask = pd.Series(False, index=pos.index)
    for lo, hi in windows:
        mask |= (pos >= lo) & (pos <= hi)
    return mask


def gtex_region(allpairs_dir: Path, tissue: str, arm: str, chrom: str,
                windows: list[tuple[int, int]]) -> pd.DataFrame:
    f = allpairs_dir / ALLPAIRS_PATTERN[arm].format(tissue=tissue, chrom=chrom)
    if not f.is_file():
        # Loud. An arm that silently produced no region would look like a null
        # colocalization result rather than a missing file, and the gate counts
        # evaluated regions -- so a typo in a filename would read as evidence.
        raise SystemExit(
            f"Missing GTEx all-pairs file for {arm} / {tissue} / {chrom}: {f}")
    # The two resources also name the phenotype column differently: eQTL
    # all-pairs use `gene_id`, sQTL all-pairs use `phenotype_id`. Read it off
    # the schema rather than assuming, and fail if neither is present.
    dataset = ds.dataset(f, format="parquet")
    names = set(dataset.schema.names)
    pheno_col = next((c for c in ("phenotype_id", "gene_id") if c in names), None)
    if pheno_col is None:
        raise SystemExit(
            f"{f} has neither phenotype_id nor gene_id: {sorted(names)}")
    # Read in record batches and keep only the window, rather than
    # materialising the chromosome and filtering afterwards. The position lives
    # inside variant_id, so there is no column for pyarrow to push a filter
    # down onto -- and an sQTL all-pairs chromosome (every intron cluster x
    # every cis variant) is large enough that the whole-table form was
    # OOM-killed. Peak memory here is one batch.
    cols = [pheno_col, "variant_id", "pval_nominal", "slope", "slope_se",
            "af", "ma_count"]
    kept: list[pd.DataFrame] = []
    for batch in dataset.to_batches(columns=cols):
        if batch.num_rows == 0:
            continue
        df = batch.to_pandas().rename(columns={pheno_col: "phenotype_id"})
        comp = split_variant(df["variant_id"])
        df = pd.concat([df.reset_index(drop=True), comp], axis=1)
        df = df[in_any_window(df["pos"], windows)]
        if df.empty:
            continue
        df = df[df["slope_se"].gt(0) & df["slope_se"].notna()]
        if not df.empty:
            kept.append(df)
    if not kept:
        return pd.DataFrame()
    out = pd.concat(kept, ignore_index=True)
    return out.rename(columns={"slope": "beta", "slope_se": "se"})


def meqtl_region(nominal_f: Path, cpgs: set[str],
                 windows: list[tuple[int, int]]) -> pd.DataFrame:
    if not nominal_f.is_file():
        return pd.DataFrame()
    tbl = ds.dataset(nominal_f, format="parquet").to_table(
        filter=pc.field("phenotype_id").isin(list(cpgs)),
        columns=["phenotype_id", "variant_id", "pval_nominal", "slope",
                 "slope_se", "af", "ma_count"])
    df = tbl.to_pandas()
    if df.empty:
        return df
    comp = split_variant(df["variant_id"])
    df = pd.concat([df.reset_index(drop=True), comp], axis=1)
    df = df[in_any_window(df["pos"], windows)]
    df = df[df["slope_se"].gt(0) & df["slope_se"].notna()]
    return df.rename(columns={"slope": "beta", "slope_se": "se"})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--chrom", required=True)
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "09_schizophrenia_risk_application" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")
    man = read_manifest(run_dir)
    cfg = load_run_config("schizophrenia", run_dir, root)
    coloc_cfg = cfg["colocalization"]
    arms = coloc_cfg["arms"]
    gtex = coloc_cfg["gtex"]
    half = int(coloc_cfg["window_bp"])

    chrom = f"chr{str(args.chrom).replace('chr', '')}"
    out_dir = run_dir / "coloc" / "regions"
    out_dir.mkdir(parents=True, exist_ok=True)
    index_rows: list[dict] = []

    loci = pd.read_csv(run_dir / "results" / "scz-loci.tsv", sep="\t")
    loci = loci[loci["chrom"].astype(str) == chrom]

    # Only loci with meQTL-supported evidence are worth a regional coloc: the
    # module's question is about risk loci that reach methylation, and running
    # every published locus in six tissues would be a genome-wide screen the
    # module does not have an FDR family for.
    support_f = run_dir / "results" / "locus-meqtl-support.tsv"
    if support_f.is_file():
        sup = pd.read_csv(support_f, sep="\t")
        if not sup.empty and "has_significant_risk_variant_cpg_meqtl" in sup:
            keep = set(sup.loc[
                sup["has_significant_risk_variant_cpg_meqtl"].astype(bool),
                "locus_id"])
            loci = loci[loci["locus_id"].isin(keep)]

    if loci.empty:
        pd.DataFrame(columns=["locus_id"]).to_csv(
            out_dir / f"{chrom}.index.tsv", sep="\t", index=False)
        print(f"{chrom}: no locus with meQTL support; nothing to prepare")
        return

    gw = pd.read_csv(run_dir / "gwas" / "gwas-sumstats-hg38.tsv.gz", sep="\t")
    gw = gw[gw["chrom"].astype(str) == chrom]
    gw = gw.rename(columns={"pos_hg38": "pos", "effect_allele": "a1",
                            "other_allele": "a2"})
    gw = add_keys(gw[["chrom", "pos", "a1", "a2", "beta", "se", "pvalue",
                      "neff", "variant_rsid", "locus_id"]])

    links = pd.read_csv(run_dir / "results" / "locus-vmr-links.tsv", sep="\t")
    burden_dir = (root / "05_cpg_meqtl_burden" / "_m" / "runs"
                  / man["upstream_cpg_meqtl_burden_run_id"])
    member = pd.read_csv(burden_dir / "results" / "tested-cpg-membership.tsv",
                         sep="\t")
    meqtl_f = (burden_dir / "results" / "meqtl" / "nominal"
               / f"{chrom}.cis_qtl_pairs.{chrom}.parquet")

    n_cases = int(cfg["locus_definition"]["gwas_n_cases"])
    n_controls = int(cfg["locus_definition"]["gwas_n_controls"])
    gwas_s = n_cases / (n_cases + n_controls)

    # coloc needs the QTL sample size, and for the meQTL arm that is this
    # cohort's donor count for this region -- taken from the accepted Module 01
    # catalog rather than hard-coded, so it cannot drift from the donors the
    # meQTL scan actually used.
    donors_f = (root / "01_vmr_catalog" / "_m" / "runs"
                / man["upstream_vmr_catalog_run_id"] / "vmr" / "donors_plink.txt")
    if not donors_f.is_file():
        raise SystemExit(f"Missing Module 01 donor list: {donors_f}")
    meqtl_n = sum(1 for line in donors_f.read_text().splitlines() if line.strip())
    if meqtl_n < 2:
        raise SystemExit(f"Implausible donor count ({meqtl_n}) in {donors_f}")

    # Loop structure matters here. Resources are the OUTER loop and loci the
    # inner one, so exactly one QTL resource is resident at a time: each is read
    # once per chromosome (re-reading all twelve tissue files per locus was the
    # dominant cost) but they are never all held together (caching all twelve
    # peaked at ~5.9 GB on chr22 and was OOM-killed). One resource, all loci,
    # then release it.
    # Loci are processed in chunks, and a resource is re-read per chunk. Peak
    # memory is set by how much of a QTL resource falls inside the windows held
    # at once, and chromosomes differ by an order of magnitude in locus count
    # (3 on chr22, 35 on chr1). Holding every window of chr1 at once would be
    # ~10x the chr22 footprint, which already peaked above 5 GB for one sQTL
    # tissue. Chunking trades a re-read per chunk for a bound that does not
    # depend on which chromosome the task drew.
    chunk_size = int(os.environ.get("SCZ_COLOC_LOCUS_CHUNK", "5"))

    n_cases = int(cfg["locus_definition"]["gwas_n_cases"])
    n_controls = int(cfg["locus_definition"]["gwas_n_controls"])
    gwas_s = n_cases / (n_cases + n_controls)

    # coloc needs the QTL sample size, and for the meQTL arm that is this
    # cohort's donor count for this region -- taken from the accepted Module 01
    # catalog rather than hard-coded, so it cannot drift from the donors the
    # meQTL scan actually used.
    donors_f = (root / "01_vmr_catalog" / "_m" / "runs"
                / man["upstream_vmr_catalog_run_id"] / "vmr" / "donors_plink.txt")
    if not donors_f.is_file():
        raise SystemExit(f"Missing Module 01 donor list: {donors_f}")
    meqtl_n = sum(1 for line in donors_f.read_text().splitlines() if line.strip())
    if meqtl_n < 2:
        raise SystemExit(f"Implausible donor count ({meqtl_n}) in {donors_f}")

    # Per-locus GWAS slices and CpG sets, computed once.
    locus_rows = [dict(r) for _, r in loci.iterrows()]
    locus_gwas: dict[str, pd.DataFrame] = {}
    locus_cpgs: dict[str, set] = {}
    for lr in locus_rows:
        lid = str(lr["locus_id"])
        lo, hi = int(lr["start"]) - half, int(lr["end"]) + half
        g_reg = gw[(gw["pos"] >= lo) & (gw["pos"] <= hi)]
        locus_gwas[lid] = g_reg[[
            "chrom", "pos", "a1", "a2", "beta", "se", "pvalue", "neff",
            "variant_rsid", "match_key", "strand_ambiguous"]].copy()
        vmr_ids = set(links.loc[links["locus_id"].astype(str) == lid, "vmr_id"])
        locus_cpgs[lid] = set(
            member.loc[member["vmr_id"].isin(vmr_ids), "cpg_id"].astype(str))

    resources: list[tuple[str, str]] = []
    for arm, dir_key in (("gtex_eqtl", "eqtl_allpairs_dir"),
                         ("gtex_sqtl", "sqtl_allpairs_dir")):
        if arms.get(arm, {}).get("enabled"):
            resources.extend((arm, t) for t in gtex["tissues"])
    if arms.get("meqtl", {}).get("enabled"):
        resources.append(("meqtl", man["region"]))

    dir_for = {"gtex_eqtl": "eqtl_allpairs_dir", "gtex_sqtl": "sqtl_allpairs_dir"}

    n_chunks = (len(locus_rows) + chunk_size - 1) // chunk_size
    for chunk_i in range(n_chunks):
        chunk = locus_rows[chunk_i * chunk_size:(chunk_i + 1) * chunk_size]
        windows = [(max(1, int(r["start"]) - half), int(r["end"]) + half)
                   for r in chunk]
        chunk_cpgs = set().union(
            *(locus_cpgs[str(r["locus_id"])] for r in chunk)) if chunk else set()
        print(f"  chunk {chunk_i + 1}/{n_chunks}: {len(chunk)} loci", flush=True)

        for arm, tissue in resources:
            if arm == "meqtl":
                if not chunk_cpgs:
                    continue
                qtl_all = meqtl_region(meqtl_f, chunk_cpgs, windows)
            else:
                qtl_all = gtex_region(Path(gtex[dir_for[arm]]), tissue, arm,
                                      chrom, windows)
            print(f"    {arm}/{tissue}: {len(qtl_all)} variants in windows",
                  flush=True)
            if qtl_all.empty:
                continue
            qtl_all = add_keys(qtl_all)

            for lr in chunk:
                locus_id = str(lr["locus_id"])
                lo, hi = int(lr["start"]) - half, int(lr["end"]) + half
                g_reg = locus_gwas[locus_id]
                if g_reg.empty:
                    continue
                q_loc = qtl_all[(qtl_all["pos"] >= lo) & (qtl_all["pos"] <= hi)]
                if arm == "meqtl":
                    if not locus_cpgs[locus_id]:
                        continue
                    q_loc = q_loc[q_loc["phenotype_id"].isin(locus_cpgs[locus_id])]
                if q_loc.empty:
                    continue

                for phenotype_id, q in q_loc.groupby("phenotype_id"):
                    m = harmonise(g_reg, q)
                    n_shared = len(m)
                    out_name = (f"{chrom}__{locus_id}__{arm}__{tissue}__"
                                f"{str(phenotype_id).replace('/', '_')}")
                    if n_shared > 0:
                        af_col = "af_qtl" if "af_qtl" in m.columns else "af"
                        keep = m[["match_key", "chrom_gwas", "pos_gwas",
                                  "a1_gwas", "a2_gwas", "beta_gwas", "se_gwas",
                                  "pvalue", "neff", "variant_rsid",
                                  "qtl_beta_aligned", "se_qtl", "pval_nominal",
                                  af_col, "allele_orientation"]].copy()
                        keep.columns = ["match_key", "chrom", "pos", "a1", "a2",
                                        "gwas_beta", "gwas_se", "gwas_pvalue",
                                        "gwas_neff", "variant_rsid", "qtl_beta",
                                        "qtl_se", "qtl_pvalue", "qtl_af",
                                        "allele_orientation"]
                        keep.to_csv(out_dir / f"{out_name}.tsv.gz", sep="\t",
                                    index=False, compression="gzip")
                    index_rows.append({
                        "region_file": f"{out_name}.tsv.gz" if n_shared else "",
                        "chrom": chrom, "locus_id": locus_id, "arm": arm,
                        "qtl_context": tissue, "phenotype_id": phenotype_id,
                        "qtl_type": "quant",
                        "n_variants_gwas": len(g_reg),
                        "n_variants_qtl": len(q),
                        "n_variants_shared": n_shared,
                        "gwas_type": cfg["locus_definition"]["gwas_type"],
                        "gwas_s": gwas_s,
                        "gwas_n": n_cases + n_controls,
                        "qtl_n": (meqtl_n if arm == "meqtl"
                                  else int(gtex["n_samples_default"])),
                        "qtl_ancestry": arms[arm]["qtl_ancestry"],
                        "ld_ancestry_matched": bool(arms[arm]["ld_ancestry_matched"]),
                        "gate_eligible": bool(arms[arm]["gate_eligible"]),
                        "arm_status": arms[arm].get("arm_status", "OK"),
                    })
            del qtl_all

    idx = pd.DataFrame(index_rows)
    idx.to_csv(out_dir / f"{chrom}.index.tsv", sep="\t", index=False)
    print(f"{chrom}: prepared {len(idx)} coloc regions across "
          f"{loci['locus_id'].nunique()} loci")


if __name__ == "__main__":
    main()
